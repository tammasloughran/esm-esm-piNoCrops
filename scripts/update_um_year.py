#!/g/data/hh5/public/apps/nci_scripts/python-analysis3
# Copyright 2020 Scott Wales
# author: Scott Wales <scott.wales@unimelb.edu.au>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import mule
import sys
import tempfile
import shutil

year = int(sys.argv[1])
restart = sys.argv[2]

mf = mule.DumpFile.from_file(restart)

mf.fixed_length_header.t1_year = year
mf.fixed_length_header.t2_year = year

for f in mf.fields:
    f.lbyr = year

mf.validate = lambda *args, **kwargs: True

temp = tempfile.NamedTemporaryFile()
mf.to_file(temp.name)

shutil.copy(temp.name, restart)
