#!/usr/bin/python3 -i

import sys
import os
import time

try:
  from gurobipy import *
  gurobi.interactive()
  _m = Model()
  del _m
except Exception as e:
  print('')
  print('%s' % e.message)
  print('')
  if e.message.find('No Gurobi license found') != -1:
    print('Running grbgetkey...')
    os.system('grbgetkey')
    print('Restart the Gurobi Interactive Shell to use a newly retrieved license file')
  time.sleep(5)
  exit(1)

version = str(gurobi.version()[0]) + '.' + \
          str(gurobi.version()[1]) + '.' + \
          str(gurobi.version()[2])
platform = gurobi.platform()

print('\nGurobi Interactive Shell (' + platform + '), Version ' + version)
print('Copyright (c) 2024, Gurobi Optimization, LLC')
print('Type "help()" for help\n')

sys.ps1 = "gurobi> "
sys.ps2 = "....... "
