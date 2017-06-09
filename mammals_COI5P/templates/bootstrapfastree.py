#!/usr/bin/env python

from subprocess import Popen, PIPE, STDOUT
import tarfile
import gzip

file = gzip.open("!{bootFile.baseName}.nw.gz", "wb")
tar = tarfile.open("!{bootFile}", "r:gz")
for tarinfo in tar:
    f=tar.extractfile(tarinfo)
    content=f.read()
    p = Popen(['FastTree','-nopr','-nosupport','-wag','-gamma'], stdout=PIPE, stdin=PIPE, stderr=PIPE)
    stdout_data = p.communicate(input=content)[0]
    file.write(stdout_data)
tar.close()
file.close()
