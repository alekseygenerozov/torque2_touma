import glob
import numpy as np
from bash_command import bash_command as bc
import shlex

import json

tags=glob.glob("*ein*")
tags=np.sort(tags)
prec_rate={}
for tag in tags:
	x=bc.bash_command("cat {0}".format(tag))
	x=float(x)
	prec_rate[tag]=x


js=json.dumps(prec_rate)
with open("prec_rate.json", "w") as outfile:
	json.dump(js, outfile)


