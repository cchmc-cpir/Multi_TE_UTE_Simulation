# Multi_TE_UTE_Simulation
### Multi-TE Lung T2* Decay Simulation

This repository contains MATLAB code for simulating multi-TE lung T2* decay using modified 3D Shepp-Logan phantoms. The simulation workflow is described below:

#### Simulation Workflow:

1. **Phantom Generation**: Modified 3D Shepp-Logan phantoms with isotropic matrix sizes of 64x64x64, 96x96x96, or 128x128x128 is created. The signal intensity pattern of the phantom resembles SNR patterns observed from murine lungs in vivo.

2. **T2* Relaxation Modeling**: T2* relaxation was incorporated into the signal decay during both the nominal echo time and the radial readout. Lung parenchyma was assigned a T2* value of 0.40 ms, and vasculature was assigned a T2* value of 3.16 ms.

3. **Signal Decay Equation**: The image signal was varied at the voxel level using the signal decay equation:

    S(k,i) = S0 * exp(-(TE(k)+(i-1)td)/T2*

    where:
    - \(S(k,I)\) is the signal at voxel \(i\) and echo time \(k\).
    - \(S0\) is the initial signal without decay.
    - \(TE\) are the nominal echo times.
    - \(i\) is the index of the number of points along the radial projections.
    - \(td\) is the dwell time.

4. **K-Space Data Generation**: For each unique combination of echo time and radial projection, k-space data were generated by performing a Fast Fourier Transform (FFT) of the signal-decayed images. Golden angle sampling was used to create radial FIDs.

5. **Noise Addition**: Gaussian noise was added separately to the real and imaginary portions of the complex k-space data to generate noisy FIDs, simulating image SNR levels typically observed in fully sampled, in-vivo mouse lung images.

6. **Image Reconstruction**: Images were reconstructed from noisy FIDs data and prescribed trajectories using Cartesian re-gridding, iterative density compensation, and FFT.

#### How to Run the Simulation:

To run the simulation, execute the file named "main" in MATLAB.

For any questions or issues, please contact Abdullah S. Bdaiwi at abdullah.bdaiwi@cchmc.org.



