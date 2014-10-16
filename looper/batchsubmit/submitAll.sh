#!/bin/bash
voms-proxy-init -voms cms -valid 240:00
condor_submit configs_V00-00-03/condor_V00-00-03_qcdmuenriched.cmd
condor_submit configs_V00-00-03/condor_V00-00-03_ttpythia.cmd
