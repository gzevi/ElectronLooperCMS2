#!/bin/bash
voms-proxy-init -voms cms -valid 240:00
condor_submit configs_V00-00-04/condor_V00-00-04_bjets.cmd
condor_submit configs_V00-00-04/condor_V00-00-04_qcdmuenriched.cmd
condor_submit configs_V00-00-04/condor_V00-00-04_ttpythia.cmd
condor_submit configs_V00-00-04/condor_V00-00-04_ttsemilep1.cmd
condor_submit configs_V00-00-04/condor_V00-00-04_ttsemilep2.cmd
