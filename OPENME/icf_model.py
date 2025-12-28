##imports
#for folder with plot images
import os
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(BASE_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)
#for math
import numpy as np
#for solving ODEs
from scipy.integrate import solve_ivp
#for plotting
import matplotlib.pyplot as plt

#for better plotting later on
plt.rcParams.update({
    "figure.dpi": 140,
    "savefig.dpi": 300,
    "font.size": 11,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
})


##constants and parameters
# Unit conversions
#Boltzmann Constant, connects temp in kelvin to energy in joules
kB = 1.380649e-23
#electron volt to joules
eV_J = 1.602176634e-19
#kiloelectron volt to joules
keV_J = 1.0e3 * eV_J
#megaelectron volt to joules
MeV_J = 1.0e6 * eV_J

#fusion energies (D-T)
#alpha particle energy in d-t reaction from MeV to Joules
Ea = 3.5 * MeV_J
#total reaction energy in d-t reaction from MeV to Joules
Er = 17.6 * MeV_J

#parameters
#initial temp in KeV
T0keV = 15.0
#time of simulation, 2 nanoseconds               
t = 2.0e-9
#radius of hot center in meters             
R = 30e-6
#hot center presumbed to be sphere, R**3 means R^3, volume calculation               
V = (4/3) * np.pi * R**3

#densities, assuming 50/50 Deuterium/Tritium
n_cm3 = 1.0e26 #presumbed total ion per cm^-3
n = n_cm3 * 1.0e6  #converts to m^-3
nD = 0.5 * n #states 50% of ions are Deuterium
nT = 0.5 * n #states 50% of ions are Tritium

#loss modeling
time = 2e-9                 #energy confinement time, smaller means energy is lost/leaks quicker. 2ns          
f_dep = 0.8                  #alpha energy amount that is deposited within hot spot center

##sigma v DT
#takes temperature and returns sigma v, sigma v being average reaction rate.
def sigma_v_DT(T_keV):
    #when temperature is <=0, then there will be 0 fusion.
    if T_keV <= 0:
        return 0.0
    #exponent of T_keV^2 allows the rate of fusion to grow when temp is higher
    #the (-19.94 / np.sqrt(T_keV)) prevents unreasonably high reaction rate at low temperatures
    return 1.1e-24 * T_keV**2 * np.exp(-19.94 / np.sqrt(T_keV))

##alpha heating power
#in watts
def alpha_heating_power(T_keV):
    #fusion reaction rate formula:
    #nD = deuteurium density
    #nT = tritium density
    #sigma_v_DT(T_keV) = reaction rate as a function of temperature
    #V = (V)olume of (h)ot (s)pot
    rate = nD * nT * sigma_v_DT(T_keV) * V
    #energy = reactions/s * joules per reaction * deposited fraction, results in Watts or J/S
    return rate * Ea * f_dep

##energy loss power
#function to compute the total energy loss power in watts
def energy_loss_power(T_keV):
    #converts keV to Kelvin
    #T_keV = avg thermal energy in keV per particle
    #T_keV * keV_J multiplies the denominator of keV on T_keV, so we get the numerator of kbT in joules
    #then we simply divide T and kbT by kb (recall that this is the Boltzmann constant), and we get T in Kelvin
    T_K = (T_keV * keV_J) / kB
    #Thermal energy of the hot spot, in joules.
    #3.0 as is common with most plasmas containing ions + electrons
    thermal_energy = 3.0 * n * kB * T_K * V
    #loss power model
    return thermal_energy / time

##derivative of temperature with respect to time, dT/dt
#derivative for ODE, t being time and y being the state vector
def dTdt(t, y):
    #extracts temperature value, as y[0] is = to temperature
    T_keV = max(y[0], 1e-6)
    #alpha heating at current temp
    #prevents negative/zero temps
    heating = alpha_heating_power(T_keV)
    #losses at current temp
    losses = energy_loss_power(T_keV)

    #conversion factor for J/s to keV/s
    heat_capacity = 3.0 * n * V * keV_J
    #temperature rises when heating is greater than lost heat (duh)
    dTdt_keV = (heating - losses) / heat_capacity
    #returning derivative
    return [dTdt_keV]

#stops the simulation when ignition marker is reached
#30 keV is the ignition marker
ign = 30.0 

##ignition event
def ignition_event(t, y):
    #triggers when at 0
    return y[0] - ign

#makes integration stop when event triggers
ignition_event.terminal = True

#triggers only when threshold is pushed
ignition_event.direction = 1


##run case
#running the model
def run_case(T0keV):
    #calls on ODE solver
    sol = solve_ivp(
        #derivative function to integrate
        dTdt,
        #time interval
        (0.0, t),
        #state vector of temperature only
        [T0keV],
        #cannot do more than 10 picosecond steps at a time, to give accurate answers.
        max_step=1e-11,
        events=ignition_event
        #if ignition threshold is reached and is exceeding, it shuts it off
        #sol.t contains times
        #sol.y[0] contains temperature
    )

    maxT = float(np.max(sol.y[0]))

    #if ignition is triggered, gives info on when
    ignited = sol.t_events[0].size > 0
    if ignited:
        t_ign = float(sol.t_events[0][0])
    else:
        t_ign = float("nan")

    return ignited, t_ign, maxT, sol

#explained with other func piece below
def main():
    ##threshold sweep
    #initial temperature values in steps of 0.5
    T0_values = np.linspace(2.0, 20.0, 37)

    #list for if ignited or not
    ignited_flags = []
    #list for ignition times in ns
    ignns = []
    #list for max temp per
    maxT_v = []

    for T0 in T0_values:
        #updates initial temperature for case being run
        ignited, t_ign, maxT, _ = run_case(T0)
        #true for ignited, else false
        ignited_flags.append(ignited)
        #takes time value in ns
        ignns.append(t_ign * 1e9 if ignited else np.nan)
        #takes max T value
        maxT_v.append(maxT)

    #list to array conversion
    #bool for if ignited or not
    ignited_flags = np.array(ignited_flags, dtype=bool)
    #float for ignition times in ns
    ignns = np.array(ignns, dtype=float)
    #float for max temp per
    maxT_v = np.array(maxT_v, dtype=float)


    if np.any(ignited_flags):
        #returns on first true value for ignited_flags
        idx_first = np.argmax(ignited_flags)
        #when ignited
        print(f"First ignition in range at T0 = {T0_values[idx_first]:.2f}keV")
    else:
        #when no ignition
        print("No ignition reached in this range.")

    ##threshold plot, plot 1
    #intial temp in steps of 0.5 x axis
    #max temp y axis
    plt.figure(figsize=(6,4))
    plt.plot(T0_values, maxT_v, marker='o', label="Max T reached")
    #arrow pointing to where ignition threshold is reached
    plt.axhline(ign, linestyle="--", label="Ignition marker (30 keV)")
    plt.xlabel("Initial Temperature (keV)")
    plt.ylabel("Maximum Temperature (keV)")
    plt.title("Ignition Threshold vs Initial Temperature (2 ns)")
    #ignition threshold @ 30keV
    if np.any(ignited_flags):
        idx_first = np.argmax(ignited_flags)
        plt.annotate("Ignition threshold",
                    xy=(T0_values[idx_first], ign),
                    xytext=(T0_values[idx_first]+1, ign-6),
                    arrowprops=dict(arrowstyle="->"))
    plt.legend()
    plt.tight_layout()
    #saving plot to folder
    plt.savefig(os.path.join(FIG_DIR, "plot1_threshold.png"), bbox_inches="tight")
    plt.show()
    plt.close()



    ##ignition speed plot, plot 2
    #intial temperature in steps of 0.5 as x axis
    #ignition time as y axis
    def plot_ignition_time(T0_values, ignns):
        plt.figure(figsize=(6,4))
        plt.plot(T0_values, ignns, marker="o", linewidth=1.5, label="t_ign")
        plt.xlabel("Initial Temperature (keV)")
        plt.ylabel("Ignition time (ns)")
        plt.title("Ignition Time vs Initial Temperature (2 ns)")
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, "plot2_ignitiontime.png"), bbox_inches="tight")
        plt.show()
        plt.close()




    ##temp records, plot 3
    #plots multiple curves of initial temperature relative to time
    #displays which initial temps ignite and at what rate relative to others
    def plot_time_records(T0_list, time_label="2 ns"):
        plt.figure(figsize=(6.6,4.2))

        for T0 in T0_list:
            ignited, t_ign, maxT, sol = run_case(T0)
            plt.plot(sol.t * 1e9, sol.y[0], linewidth=1.6, label=f"{T0:.1f} keV")

        plt.axhline(ign, linestyle="--", linewidth=1.2, label="Ignition marker (30 keV)")
        plt.xlabel("Time (ns)")
        plt.ylabel("Temperature (keV)")
        plt.title(f"Temperature Records ({time_label})")
        plt.legend(ncol=2)
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, "plot3_temprecords.png"), bbox_inches="tight")
        plt.show()
        plt.close()



    #call plot 2
    plot_ignition_time(T0_values, ignns)
    #call plot 3
    plot_time_records([5.0, 10.0, 10.5, 11.0, 15.0], time_label="2 ns")

    print(f"energy confinement time = {time*1e9:.2f}ns, R = {R*1e6:.1f}µm, n = {n:.2e}m^-3")
    if np.any(ignited_flags):
        print(f"First ignition in sweep: T0 ≈ {T0_values[idx_first]:.2f}keV, t_ign ≈ {ignns[idx_first]:.3f}ns")

#allows the data to not run sweeps and plots if imported to other script later
#irrelevant for this project
if __name__ == "__main__":
    main()