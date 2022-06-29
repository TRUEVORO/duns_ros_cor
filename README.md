# Flow of gas and liquid mixtures in wells #

This code can help you calculate pressure and temperature in pipeline

#### Current Version: 0.0.01

### General structure of the project

* lib - modules for calculation process
* output_data - output data files

## Example
'''python
# Input data
p_wh = 20 * 101325
t_wh = 293
md_vdp = 1800
d_tub = 0.1524
roughness = 18.288 * 1e-6
theta_deg = 90
q_fluid_m3s = np.linspace(10, 400, 40)
gor = 100
rho_lrc_kgm3 = 762.64
rho_grc_kgm3 = 94.19
mul_rc_kgm3 = 0.97 * 1e-3
mug_rc_kgm3 = 0.016 * 1e-3
sigma_l_nm = 8.41 * 1e-3
temp_grad = 100

# Class initialization
ppl = pipe.Pipeline()

# VLP results
q_liq_m3s, p_wfs = ppl.calc_vlp(
  p_wh=p_wh,
  t_wh=t_wh,
  h0=0,
  md_vdp=md_vdp,
  d_tub=d_tub,
  roughness=roughness,
  theta_deg=theta_deg,
  q_fluid_m3s=q_fluid_m3s / 86400,
  gor=gor,
  rho_lrc_kgm3=rho_lrc_kgm3,
  rho_grc_kgm3=rho_grc_kgm3,
  mul_rc_kgm3=mul_rc_kgm3,
  mug_rc_kgm3=mug_rc_kgm3,
  sigma_l_nm=sigma_l_nm,
  temp_grad=temp_grad
  )
q_liq = q_liq_m3s * 86400

# Saving output data
with open('output_data/output.json', 'w') as ouput_file:
  json.dump({"q_liq": q_liq.tolist(), "p_wf": p_wfs.tolist()}, ouput_file)
  
# Graphs
plt.plot(q_liq, p_wfs)
plt.title('VLP')
plt.xlabel('Q_liq, m3/day')
plt.ylabel('P_wf, Pa')
plt.show()
'''
