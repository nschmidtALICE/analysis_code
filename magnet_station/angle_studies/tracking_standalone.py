import numpy as np
import ROOT
from math import fabs

# Define the magnetic field strength is gaussian along the z-axis, with at peak at 450cm and points in the y-direction
def magnetic_field(x, y, z):
    B0 = 1.1 #Tesla
    sigma = 1000.0 #to be determined
    Bx = 0.0
    By = B0 * np.exp(-(z-4500)**2 / (2.0 * sigma**2))
    # print('By', By)
    Bz = 0
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
    print('state[2]', state[2], 'target_z', target_z, 'step_size', step_size, 'target_z - state[2]', target_z - state[2])
    while (target_z - state[2]) > step_size and state[2] > 0:
    # while abs(state[2] - target_z) > step_size:
        # print('state', state)
        # print('state[2] - target_z', state[2] - target_z)
        B = magnetic_field(state[0], state[1], state[2])
        k1 = compute_derivatives(state, charge, B)
        k2 = compute_derivatives(state + 0.5 * step_size * k1, charge, B)
        k3 = compute_derivatives(state + 0.5 * step_size * k2, charge, B)
        k4 = compute_derivatives(state + step_size * k3, charge, B)
        state += (step_size / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    return state

# Example usage
def main():
    # Example input file
    inputfile = '/home/niviths/Downloads/magnetStationSims/20250214_backup_all/minimumBias_MS_MagDown_10015plus.root'
    fin = ROOT.TFile(inputfile)
    ntup = fin.Get('ntup')

    hmatch = ROOT.TH1F('hmatch', '', 100, 0, 200)

    hXZextrapolation = ROOT.TH2F('hXZextrapolation', '', 100, 3500, 8000, 100, -2500, 2500)

    numevt = 0
    #print number of events in ntup
    print(ntup.GetEntries())
    for evt in ntup:
        numevt += 1
        if numevt > 4:
            break
        print('evt', numevt)
        for i, ut_vz in enumerate(evt.ut_vz):
            print('track number', i)
            if evt.p[i] > 5000:
                print('skipping track with p > 5000')
                continue
            charge = evt.pid[i] / abs(evt.pid[i])
            track_v0 = np.array([evt.ut_vx[i], evt.ut_vy[i], ut_vz,
                                 evt.ut_tx[i], evt.ut_ty[i], evt.p[i]])
                                #  evt.ut_tx[i] * evt.p[i], evt.ut_ty[i] * evt.p[i], evt.p[i]])
            z_old = ut_vz
            track_v1 = np.array([evt.ut_vx[i], evt.ut_vy[i], ut_vz,
                                 evt.ut_tx[i] * 0.9, evt.ut_ty[i] * 0.9, evt.p[i] * 0.9])
            track_v2 = np.array([evt.ut_vx[i], evt.ut_vy[i], ut_vz,
                                 evt.ut_tx[i] * 1.1, evt.ut_ty[i] * 1.1, evt.p[i] * 1.1])

            z_finals = []
            y_finals = []
            x_finals = []
            ntrack = 0
            z_target = 0.0
            for track_v in [track_v1, track_v2]:
                ntrack += 1
                # print('track', ntrack)
                loopnum = 0
                while (fabs(track_v[0]) < 2500.0) and (fabs(track_v[1]) < 1500.0) and (z_old < 7400.0) and fabs(track_v[3] < 0.9):
                    z_target = z_old + 10.0
                    # print(z_target)
                    loopnum += 1
                    # print('loopnum', loopnum)
                    track_v = propagate_track_rk(track_v, charge, 1.0, z_target)
                    z_old = z_target
                # if z_target < 7000:
                #     print('z', z_target, 'x', track_v[0], 'z', track_v[2])
                z_finals.append(z_target)
                y_finals.append(track_v[1])
                x_finals.append(track_v[0])

            #print z finals and y finals as well as the initial track state
            # print('track_v0', np.ceil(track_v0))
            # print('evt.ut_tx[i], evt.ut_ty[i]', evt.ut_tx[i], evt.ut_ty[i])
            print('x_finals', np.ceil(x_finals))
            # print('y_finals', y_finals)
            print('z_finals', np.ceil(z_finals))

            #use x_finals and z_finals to fill the histogram
            for x, z in zip(x_finals, z_finals):
                print('x', x, 'z', z)
                hXZextrapolation.Fill(x, z)


            nhits_plane1 = 0
            for ims, ms_vz in enumerate(evt.ms_vz):
                if (ms_vz > z_finals[0]) & (ms_vz < z_finals[1]) & (evt.ms_vy[ims] > min(y_finals)) & (evt.ms_vy[ims] < max(y_finals)):
                    nhits_plane1 += 1

            hmatch.Fill(nhits_plane1 / 4)

    #plot the histogram
    c1 = ROOT.TCanvas('c1', 'c1', 800, 600)
    hXZextrapolation.Draw('colz')
    c1.SaveAs('extrapolation.png')

    fhisto = ROOT.TFile('histo_matches.root', 'recreate')
    hmatch.Write()
    hXZextrapolation.Write()
    fhisto.Close()

if __name__ == "__main__":
    main()