# run using the command
# lb-run Moore/latest python fast_tracking.py

import ROOT
import numpy as np
from math import fabs
import os
import GaudiPython as GP
from GaudiPython.Bindings import gbl, AppMgr, iAlgTool
from Configurables import LHCbApp
# from Configurables import ToolSvc 
ROOT.gInterpreter.Declare('#include <GaudiKernel/SystemOfUnits.h>')
Units = gbl.Gaudi.Units

#name = os.environ['NAME']
name = '/eos/user/n/nschmidt/MagnetStation/minimumBias_MS_MagDown_1900plus.root'

app = LHCbApp(
    DataType="Upgrade",
    EvtMax=10,
    Simulation=True,
    DDDBtag   = "dddb-20231017",
    CondDBtag = "sim-20231017-vc-md100")

# Upgrade DBs
# CondDB().Upgrade = True
appMgr = GP.AppMgr()
appMgr.initialize()
toolSvc = appMgr.toolSvc()
LHCb = gbl.LHCb
m_extrapolator = toolSvc.create('TrackRungeKuttaExtrapolator',
                                interface='ITrackExtrapolator')

hmatch = ROOT.TH1F('hmatch','',100, 0, 200)

fin = ROOT.TFile(name)
ntup = fin.Get('ntup')
zbars = np.linspace(3500, 7000, 700)
ybars = np.array([50, 80, 140, 230, 350, 500, 680, 970, 1200])

numevt = 0

dummy_matrix = ROOT.Math.SMatrixSym5D()

# Define the magnetic field (example: constant field in the z direction)
def magnetic_field(x, y, z):
    Bx = 0.0
    By = 0.0
    Bz = 1.0  # Tesla
    return np.array([Bx, By, Bz])

# Function to compute the derivatives of the track state
def compute_derivatives(state, charge, B):
    x, y, z, px, py, pz = state
    momentum = np.sqrt(px**2 + py**2 + pz**2)
    vx = px / momentum
    vy = py / momentum
    vz = pz / momentum
    ax = charge * (vy * B[2] - vz * B[1])
    ay = charge * (vz * B[0] - vx * B[2])
    az = charge * (vx * B[1] - vy * B[0])
    return np.array([vx, vy, vz, ax, ay, az])

# Runge-Kutta 4th order method to propagate the track state
def propagate_track_rk(state, charge, step_size, target_z):
    while abs(state[2] - target_z) > step_size:
        B = magnetic_field(state[0], state[1], state[2])
        k1 = compute_derivatives(state, charge, B)
        k2 = compute_derivatives(state + 0.5 * step_size * k1, charge, B)
        k3 = compute_derivatives(state + 0.5 * step_size * k2, charge, B)
        k4 = compute_derivatives(state + step_size * k3, charge, B)
        state += (step_size / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    return state

for evt in ntup:
    numevt += 1
    #print event number
    print('evt', numevt)
    for i, ut_vz in enumerate(evt.ut_vz):
        charge = evt.pid[i] / abs(evt.pid[i])
        track_v0 = np.array([evt.ut_vx[i] * Units.mm, evt.ut_vy[i] * Units.mm, ut_vz * Units.mm,
                             evt.ut_tx[i] * evt.p[i], evt.ut_ty[i] * evt.p[i], evt.p[i]])
        z_old = ut_vz
        track_v1 = np.array([evt.ut_vx[i] * Units.mm, evt.ut_vy[i] * Units.mm, ut_vz * Units.mm,
                             evt.ut_tx[i] * evt.p[i] * 0.9, evt.ut_ty[i] * evt.p[i] * 0.9, evt.p[i] * 0.9])
        track_v2 = np.array([evt.ut_vx[i] * Units.mm, evt.ut_vy[i] * Units.mm, ut_vz * Units.mm,
                             evt.ut_tx[i] * evt.p[i] * 1.1, evt.ut_ty[i] * evt.p[i] * 1.1, evt.p[i] * 1.1])

        z_finals = []
        y_finals = []
        ntrack = 0
        for track_v in [track_v1, track_v2]:
            ntrack += 1
            print('track', ntrack)
            while (fabs(track_v[0]) < 2500.0) and (fabs(track_v[1]) < 1500.0) and (z_old < 7400.0) and (track_v[2] < 1):
                z_target = z_old + 10.0
                print(z_target)
                track_v = propagate_track_rk(track_v, charge, 1.0, z_target)
                z_old = z_target
            z_finals.append(z_target)
            y_finals.append(track_v[1])

        nhits_plane1 = 0
        for ims, ms_vz in enumerate(evt.ms_vz):
            if (ms_vz > z_finals[0]) & (ms_vz < z_finals[1]) & (evt.ms_vy[ims] > min(y_finals)) & (evt.ms_vy[ims] < max(y_finals)):
                nhits_plane1 += 1

        hmatch.Fill(nhits_plane1 / 4)

fhisto = ROOT.TFile('histo_matches.root', 'recreate')
hmatch.Write()
fhisto.Close()