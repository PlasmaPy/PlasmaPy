from dispersionsolver import two_fluid_dispersion_solution as tfds

k = np.logspace(2, 7, 1000)
zeta, v_A, beta, c_s = tfds(n=5, B=5.e-5, T_i=1.e5, theta=np.pi/4, k=k)

omega_fast     = np.sqrt(zeta['fast_mode'] * k**2 * v_A**2)
omega_alfven   = np.sqrt(zeta['alfven_mode'] * k**2 * v_A**2)
omega_acoustic = np.sqrt(zeta['acoustic_mode'] * k**2 * v_A**2)

plt.clf()

plt.plot(k, omega_fast, 'r.', ms=5, label='fast')
plt.plot(k, omega_alfven, 'b.', ms=5, label='alfven')
plt.plot(k, omega_acoustic, 'g.', ms=5, label='acoustic')
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.show()

print('omega_fast=%s,%s' %(np.nanmin(omega_fast), np.nanmax(omega_fast)))
print('omega_alfven=%s,%s' %(np.nanmin(omega_alfven), np.nanmax(omega_alfven)))
print('omega_acoustic=%s,%s' %(np.nanmin(omega_acoustic), np.nanmax(omega_acoustic)))

