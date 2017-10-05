import os
import sys
import re
import gzip
import yaml
import json
import argparse


class FastqFilename:
    pattern = re.compile(r"(?P<name>.*)(?P<extension>\.fastq|\.fastq\.gz|\.fq|fq\.gz)$")
    match = pattern.match


class IlluminaFastqFilename:
    pattern = re.compile(r"(?P<name>.+?)_?(?P<barcode>[NACTG]{3,30})?_?(?P<lane>L\d{3})?(_)?(?P<read>R\d)?_?(?P<set>\d{3})(?P<extension>\.fastq|\.fastq\.gz)$")
    match = pattern.match


class IlluminaSeqidV1:
    pattern = re.compile(r"@(?P<instrument>[a-zA-Z0-9_-]*):(?P<lane>\d*):(?P<tile>\d*):(?P<x_pos>\d*):(?P<y_pos>\d*)(?P<index_number>#\d|[NACTG]*)\/(?P<read>\d)")
    match = pattern.match


class IlluminaSeqidV2:
    pattern = re.compile(r"@(?P<instrument>[a-zA-Z0-9_-]*):(?P<run_number>\d*):(?P<flowcellID>[a-zA-Z0-9]*):(?P<lane>\d*):(?P<tile>\d*):(?P<x_pos>\d*):(?P<y_pos>\d*)\s(?P<read>\d*):(?P<is_filtered>[YN]):(?P<control_number>\d*):(?P<index_sequence>[NACTG]*)")
    match = pattern.match


filename_patterns = (IlluminaFastqFilename, FastqFilename)
seqid_patterns = (IlluminaSeqidV2, IlluminaSeqidV1)


class Fastq:
    def __init__(self, path):
        """ Fastq container object that resolves metadata from the given path.
        Attempts to extract Name, Lane, and Read from filename, and FCID from
        the first line of the FASTQ.
        :param path: Path to a FASTQ file.
        """
        self.path = path
        self.filename = os.path.basename(self.path)
        if not FastqFilename.match(self.filename):
            raise ValueError(' does not match fastq pattern')

        # Parse the filename
        for pattern in filename_patterns:
            groups = pattern.match(self.filename)
            if groups:
                gd = groups.groupdict()
                self.filename_pattern = pattern.__name__
                self.name = gd['name']
                self.lane = gd.get('lane') or '1'
                self.read = gd.get('read') or '1'
                self.read = int(self.read.strip('R'))

        # Read the first seqid
        if self.path.endswith('gz'):
            self.gzipped = True
            with gzip.open(self.path, 'rb') as fq:
                self.seqid = fq.readline().decode().strip()
        else:
            self.gzipped = False
            with open(self.path, 'r') as fq:
                self.seqid = fq.readline().strip()

        # Parse the seqid
        for pattern in seqid_patterns:
            groups = pattern.match(self.seqid)
            if groups:
                gd = groups.groupdict()
                self.seqid_pattern = pattern.__name__
                self.instrument = gd.get('instrument')
                self.run_number = gd.get('run_number')
                self.fcid = gd.get('flowcellID')
                self.lane = gd.get('lane') or self.lane
        else:
            self.seqid_pattern = None
            self.fcid = 'Unknown'
        # Beginnings of a read group tag for this fastq
        self.readgroup = "{}{}".format(self.fcid, self.lane)

    def __str__(self):
        return yaml.dump(self.__dict__, default_flow_style=False)

    def length(self, openfn=None, mode=None):
        """ Get the total number of reads"""
        if openfn is None or mode is None:
            if self.gzipped:
                openfn = gzip.open
                mode = 'rb'
            else:
                openfn = open
                mode = 'r'

        i = 0
        with openfn(self.path, mode) as fq:
            for i, l in enumerate(fq):
                pass
        return (i + 1)/4


def as_yaml(mapping):
    """ Returns mapping object as a pretty yaml string"""
    return yaml.dump(mapping, default_flow_style=False)


def as_json(mapping):
    """ Returns mapping object as a pretty json string"""
    return json.dumps(mapping, indent=4)


def load_json(path):
    with open(path, 'r') as fp:
        return json.load(fp)


def load_yaml(path):
    with open(path, 'r') as fp:
        return yaml.load(fp.read())


def filter_paths(paths, filter):
    results = [path for path in paths if filter.match(path)]
    return results


def scan_dir(path):
    """ Scan a directory and add any fastqs to samples """
    fastq_paths = []
    abspath = os.path.abspath(os.path.expanduser(path))
    for dp, dn, filenames in os.walk(abspath):
        for fn in filenames:
            if FastqFilename.match(fn):
                p = os.path.join(dp, fn)
                fastq_paths.append(p)
    return fastq_paths


def build_samples(fastqs):
    samples = {}
    for fastq in fastqs:
        if fastq.name not in samples:
            samples[fastq.name] = {'name': fastq.name}
        if 'name' not in samples[fastq.name]:
            samples[fastq.name]['name'] = fastq.name
        if 'readgroups' not in samples[fastq.name]:
            samples[fastq.name]['readgroups'] = {}
        if fastq.readgroup not in samples[fastq.name]['readgroups']:
            samples[fastq.name]['readgroups'][fastq.readgroup] = [None]
        while len(samples[fastq.name]['readgroups'][fastq.readgroup]) < fastq.read:
            samples[fastq.name]['readgroups'][fastq.readgroup].append(None)
        # Making an assumption here that no fastq will ever be labeled R0
        # This is probably a bad assumption, but I will need to write a container
        # for the reads that can be 0 or 1 base indexed depending on what is added?
        samples[fastq.name]['readgroups'][fastq.readgroup][fastq.read - 1] = fastq.path
    return samples


def print_samplesheet(samplesheet, fp, format='json'):
    if format == 'json':
        serializer = as_json
    elif format == 'yaml':
        serializer = as_yaml
    else:
        raise ValueError('Cant print samplesheets with "%s"' % format)
    fp.write(serializer(samplesheet))
    fp.flush()
    fp.close()


def load_samplesheet(path, format='json'):
    if format == 'yaml':
        loader = load_yaml
    elif format == 'json':
        loader = load_json
    else:
        raise ValueError('Cant load samplesheets with "%s"' % format)

    samplesheet = loader(path)
    return samplesheet


def get_args_discover(args=None):
    """ Construct the argument parser """
    parser = argparse.ArgumentParser(description="""
       Discover FASTQ files and create a samplesheet.
       Searches path(s) for FASTQ files and extracts sample
       information from the filename and sequence identifiers to
       construct a samplesheet. The resulting samplesheet is printed
       as YAML or JSON document to stdout.
       """)

    parser.add_argument('path', nargs='+',
                        help='Path to search')
    parser.add_argument('-a', '--append', default=None,
                        help='Append samples to existing samplesheet')
    parser.add_argument('-l', '--loader', default='yaml',
                        help='Format of existing samplesheet')
    parser.add_argument('-f', '--filter', default=None,
                        help='Regex pattern to filter matching files')
    parser.add_argument('-o', '--out', default=None, type=argparse.FileType('w'),
                        help='Output path, "-" for stdout')
    parser.add_argument('--format', default='yaml',
                        help='Output format')
    parser.add_argument('-q', '--quiet', action='store_true', default=False,
                        help='Suppress the logging on stderr')

    args = parser.parse_args(args)

    # Convert the filter to a regex object
    if args.filter:
        pat = re.compile(args.filter)
        args.filter = pat

    return args


def discover(args=None):
    args = get_args_discover(args)

    # Search for fastqs in each path given
    fastq_paths = []
    for path in args.path:
        fastq_paths += scan_dir(path)

    # Filter against a pattern
    if args.filter:
        fastq_paths = filter_paths(fastq_paths, args.filter)

    # Make Fastq objects
    fastqs = [Fastq(path) for path in fastq_paths]

    # Log a summary of what we found
    if not args.quiet:
        sys.stderr.write('Filename\tSequenceID\tPath\n')
        for fastq in fastqs:
            msg = '{filename}\t{seqid}\t{path}\n'.format(
                filename=fastq.filename_pattern,
                seqid=fastq.seqid_pattern,
                path=fastq.path
            )
            sys.stderr.write(msg)
        sys.stderr.flush()

    # Print a samplesheet if we want to
    if args.out is not None:
        samples = build_samples(fastqs)
        if args.append:
            samplesheet = load_samplesheet(args.append, args.loader)
            if 'samples' not in samplesheet:
                samplesheet['samples'] = samples
        else:
            samplesheet = {'samples': samples}

        print_samplesheet(samplesheet, args.out, args.format)


def get_args_check(args=None):
    parser = argparse.ArgumentParser(description="Investigate fastq files")
    parser.add_argument('path', nargs='+')
    args = parser.parse_args(args)
    return args


def check(args=None):
    args = get_args_check(args)
    for path in args.path:
        fastq = Fastq(path)
        print(fastq)


def main(args=None):
    parser = argparse.ArgumentParser(
        description="Fastq support for tidy sample sheets",
        add_help=False
    )
    parser.add_argument('action', choices=['discover', 'check'])
    known, args = parser.parse_known_args(args)

    action = known.action

    if action == 'discover':
        discover(args=args)
    elif action == 'check':
        check(args=args)


if __name__ == '__main__':
    main()
