#!/usr/bin/env python

from subprocess import Popen, PIPE, STDOUT
import tarfile
import gzip
import sys

file = gzip.open("boot.nw.gz", "wb")
tar = tarfile.open("boot.tar", "r")
for tarinfo in tar:
    f=tar.extractfile(tarinfo)
    content=gzip.GzipFile(fileobj=f).read()
    p = Popen(['FastTree','-nopr','-nosupport','-wag','-gamma'], stdout=PIPE, stdin=PIPE, stderr=PIPE)
    stdout_data , err = p.communicate(input=content)
    file.write(stdout_data)
    if err:
        sys.stderr.write(err)
tar.close()
file.close()
