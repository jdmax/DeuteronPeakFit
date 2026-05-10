# Example of deuteron fit usage. Plots an example event from Run Group C.

from deuteron_fit import fit
import matplotlib.pyplot as plt
import json

with open("example_data/example_data.json", "r") as event:
    for line in event:
        json_dict = json.loads(line.rstrip('\n|\r'))
        signal = json_dict['fitsub']
        freqs = json_dict['freq_list']

initial_params = {
    'A':   0.03,
    'G':  -0.00003,
    'r':   1.2,      # r > 1 is positive polarization
    'wQ':  0.027,
    'wL':  32.69,
    'eta': -0.02,
    'xi':  -0.001,
}

result = fit(freqs, signal, initial_params)

r      = result.params['r'].value
r_err  = result.params['r'].stderr
denom  = r**2 + r + 1

P   = (r**2 - 1) / denom
Pzz = (r - 1)**2 / denom

if r_err is not None:
    P_err   = (r**2 + 4*r + 1) / denom**2 * r_err
    Pzz_err = 3 * abs(r**2 - 1)  / denom**2 * r_err
    print(f"Vector polarization:  P   = {P*100:.2f} ± {P_err*100:.2f}%")
    print(f"Tensor polarization:  Pzz = {Pzz*100:.2f} ± {Pzz_err*100:.2f}%")
else:
    print(f"Vector polarization:  P   = {P*100:.2f}% (uncertainty unavailable)")
    print(f"Tensor polarization:  Pzz = {Pzz*100:.2f}% (uncertainty unavailable)")
    print("Warning: lmfit could not compute parameter uncertainties.")

plt.plot(freqs, signal, label='Signal')
plt.plot(freqs, result.best_fit, '-r', label='Fit')
plt.grid()
plt.legend()
plt.title(f'NMR Signal  |  P = {P*100:.2f}%  Pzz = {Pzz*100:.2f}%')
plt.show()
