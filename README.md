0D ICF Model



I. Overview


This project presents a zero-dimensional model of an inertial confinement fusion (or, ICF) hot spot. The model tracks the evolution of the hot spot and focuses on ignition behavior under fixed conditions, such as shape and density. The primary objective is to identify the initial temperature threshold for ignition and to characterize how ignition time depends on initial conditions.

The hot spot is presumed to be a spherical plasma composed of a 50/50 deuterium/tritium mixture. The state variable evolved in time is the thermal energy per particle, expressed as kT in units of keV (kiloelectron volts). This convention is standard in plasma physics and allows fusion reactivity and heating terms to be expressed directly as functions of temperature.

Alpha particle heating from D–T fusion reactions is included. Energy losses are modeled on a scale of time and represent all different loss mechanisms. An ignition event is used to terminate the simulation once the hot spot temperature reaches 30 keV.


II. Results and Discussion


Figure 1 illustrates the maximum temperature achieved during the simulation as a function of the initial temperature for a confinement time of 2ns. The model displays the ignition threshold at an initial temperature of approximately 10–11 keV. Below this threshold, energy is lost quicker than it can be gained, and cooling occurs. Above this threshold, alpha-particle self-heating exceeds heat losses, causing a thermal runaway.

Figure 2 shows the ignition time as a function of initial temperature for cases that successfully ignite. As the initial temperature becomes larger, the ignition time becomes smaller. The decrease in time is non-linear, approaching 0 ignition time exponentially as the initial temperature becomes larger.

Figure 3 presents representative temperature responses for several initial temperatures below, near, and above the ignition threshold. Sub-threshold cases show the plasma cooling, at a slower rate the nearer to the threshold they are. Cases well above the threshold demonstrate rapid energy runaway due to alpha-particle self-heating. These trajectories illustrate the difference in behavior at different proximity and position to the threshold.


III. Explicit Assumptions

The model is based on the following assumptions:

-The hot spot is spatially uniform (0D)

-Density and volume remain constant

-Only deuterium–tritium fusion reactions are considered

-Alpha particles provide the sole source of fusion self-heating

-Energy losses are modeled using a single confinement time of 2ns

These assumptions simplify the physical system while retaining the essential nonlinear feedback responsible for ignition.


IV. Summary


Despite its simplicity, the model captures the prominent features of ICF ignition, including a sharp ignition threshold and strong sensitivity of ignition time relative to initial temperature. The results demonstrate how alpha-particle self-heating can overcome energy losses above a certain initial temperature, leading to thermal runaway and ignition under fixed confinement conditions.

V. Running the program

Read file "requirements.txt". Open code, and run. Plots will be saved to "figures" folder.
