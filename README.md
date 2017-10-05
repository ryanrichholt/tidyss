# tidyss
A collection of tools for making sample sheets better

To install:

```shell
pip install --user --upgrade git+git://github.com/ryanrichholt/tidyss.git
```

The only module available right now is fastq.py described below

## Fastqs

You can easily parse attributes from fastq files with `tidyss-fastq check`:

```shell
$ tidyss-fastq check mysample.fastq

barcode: NGATTACA
control_number: '0'
fcid: HXXXXXXXX
filename: mysample.fastq
filename_pattern: FastqFilename
gzipped: false
instrument: K01234
is_filtered: N
lane: '6'
name: mysample
path: /path/to/my.fastq
read: null
readgroup: HLY7YBBXX.6
run_number: '82'
seqid: '@K01234:82:HXXXXXXXX:6:1101:2027:1068 1:N:0:NGATTACA'
seqid_pattern: IlluminaSeqIdV2
```

You can search for fastqs and see some attributes with `tidyss-fastq discover`:

```shell
$ tidyss-fastq discover example/
Filename	SequenceID	Path
IlluminaFastqFilename	IlluminaSeqIdV2	/scratch/example/sample1_R1.fastq.gz
IlluminaFastqFilename	IlluminaSeqIdV2	/scratch/example/sample1_R2.fastq.gz
```

And you can also create tidy sample sheets `tidyss-fastq discover -o - example/`

```shell
$ tidyss-fastq discover example/ -o -
```
```yaml
samples:
  sample1:
    name: sample1
    readgroups:
      HXXXXXXX.2:
        '1': /scratch/example/sample1.fastq.gz
  sample1:
    name: sample1
    readgroups:
      HXXXXXXX.2:
        '2': /scratch/example/sample1.fastq.gz
```

