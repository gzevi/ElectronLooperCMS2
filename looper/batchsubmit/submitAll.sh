#!/bin/bash
voms-proxy-init -voms cms -valid 240:00
condor_submit configs_V00-00-06/condor_V00-00-06_qcdmuenriched.cmd
condor_submit configs_V00-00-06/condor_V00-00-06_ttsemilep1.cmd
condor_submit configs_V00-00-06/condor_V00-00-06_ttsemilep2.cmd
condor_submit configs_V00-00-06/test.cmd
