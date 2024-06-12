
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]

## Modeling the transmission of respiratory diseases in hospital setting

### Context 

Respiratory infections are a major public health issue in hospital settings. They can be transmitted from patient to patient, as well as between patients and healthcare workers, resulting in nosocomial or hospital-acquired infections. These infections contribute to the worsening of patients' health and complicate their care.

Pathogens can be transmitted in different ways:
- Large droplets emitted during speech or coughing, which can lead to short-distance transmission
- Small aerosols emitted during breathing, which can remain suspended in the air and cause long-distance transmission

A thorough understanding of the respective roles of these transmission modes is essential for developing effective strategies to prevent nosocomial respiratory infections. To date, no model simultaneously integrates these two modes of transmission, and the relative importance of each remains to be evaluated.

The aim of this project is to develop a mathematical model of respiratory pathogen transmission in a hospital ward, combining inter-individual transmission (via an epidemiological model) and airborne transmission (via a biophysical model).


## Instructions

### Dependencies
<details>
  <summary> R </summary>

* dplyr (1.1.4)
* tidyr (1.3.1)
* tidyverse (2.0.0)
* lubridate (1.9.3)
* ggplot2 (3.5.0)
* ggh4x (0.2.8)
* hrbrthemes (0.8.7)
* viridis (0.6.5)
* igraph (2.0.3)
* Rcpp (1.0.12)
* epicontacts (1.1.4)
* gganimate (1.0.9)
* gifski (1.12.0)
</details>

<details>
  <summary> c++ </summary>

  * c++ 11

  </details>


### Nods-Cov-2 scripts (R/nodscov2/)
The project is divided into 4 dependent R scripts:

<details>
  <summary> 1. interaction-nodscov2.Rmd </summary>

  This script loads Nods-Cov-2 local data (private), filters it to keep only intensive care units, then does an analysis of the distribution of individuals' categories, type of interactions etc.
  this script led us to choose (for now) Raymond Poincaré hospital's intensive care unit (Garches, France).
  
  An analysis on healthcare worker cumulative time spent interacting with patients is also done.

</details>

<details>
  <summary> 2. localization-nodscov2.Rmd </summary>
  
  This script loads Nods-Cov-2 local data (private), filters it to keep only data relating to Raymond Poincaré hospital's reanimation ward and reconstructs the localizations of individuals using interactions, inferred healthcare workers' shifts, individual category and comportemental rules.
  Localizations and other important data are saved as _/out/loc-nodscov2/localization-nodscov2.RData_.
  
  An analysis on healthcare worker cumulative time spent in the corridor and HCW's restroom is also done.

  Details on the assumptions made in a future document.

  </details>


<details>
  <summary> 3. model-nodscov2.Rmd </summary>
  This script loads Nods-Cov-2 localization data (& other related objects) from _/out/loc-nodscov2/localization-nodscov2.RData_, defines the objects needed for the simulation (parameters etc...) then compiles model-nodscov2.cpp (and its dependency model-nodscov2_fun.cpp) that is used for the epidemic simulation.
  <details>
    <summary> model-nodscov2_fun.cpp breakdown </summary>
    This c++ file contains parameters and functions used by model-nodscov2.cpp, here is a non-exhaustive list:
    <details>
      <summary>Update_environment</summary>
      This function update the viral load in each room for each time step. It is taking into account viral inactivation of the previous time step viral load and different scenarios depending on the hospital rooms (FUTURE IMPLEMENTATION).
      For the time step t, the steps are the following: 
      1. We apply an exponential decay $exp(-\mu \delta t)$ to the viral load at t-1 in each room $$ \forall t > 2, \, E_{k}(t) = \mathbf{E_{k}(t-1) e^{- \mu \delta t}} + \nu \delta t \sum_{j \, \epsilon \, I} {1_{\left \{ S(j, k, t_{i - 1})= 1 \right \}}} $$
    </details>
    <details>
      <summary>Update_environment</summary>
      This function update the viral load in each room for each time step. It is taking into account viral inactivation of the previous time step viral load and different scenarios depending on the hospital rooms (FUTURE IMPLEMENTATION).
      For the time step t, the steps are the following: 
      1. We apply an exponential decay $exp(-\mu \delta t)$ to the viral load at t-1 in each room $$ \forall t > 2, \, E_{k}(t) = \mathbf{E_{k}(t-1) e^{- \mu \delta t}} + \nu \delta t \sum_{j \, \epsilon \, I} {1_{\left \{ S(j, k, t_{i - 1})= 1 \right \}}} $$
    </details>


  </details>
  Simulation results, localizations, interactions and parameters used are sad as _/out/sim-nodscov2/&lt;id_sim&gt;-simulation-nodscov2.RData_
</details>

<details>
  <summary> 4. visualization-nodscov2.Rmd </summary>
  
  This script loads Nods-Cov-2 epidemic simulation results (& other related objects such as localizations, simulations, number of days simulated and model parameters used) from _/out/sim-nodscov2/&lt;id_sim&gt;-simulation-nodscov2.RData_ then proceeds to generate multiple plots such as:

  * SEIR population dynamic over time
  * Viral load in each room over time
  * Lambda_e (environmental) and Lambda_c (close-contact) for patients (PA) and healthcare workers (PE) over time
  * Force of Infection (FOI) for patients (PA) and healthcare workers (PE) over time
  * Infection network (thanks to the <a href='https://cran.r-project.org/web/packages/epicontacts/index.html'>epicontacts</a> R package)
  * Individual's trajectories over time

</details>



## License

Distributed under the GNU General Public License. See `COPYING` for more information.


## Contact

Project Link: [https://github.com/nollive/projet-pacri](https://github.com/nollive/projet-pacri)


<!-- MARKDOWN LINKS & IMAGES -->
[contributors-shield]: https://img.shields.io/github/contributors/nollive/projet-pacri.svg?style=for-the-badge
[contributors-url]: https://github.com/nollive/projet-pacri/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/nollive/projet-pacri.svg?style=for-the-badge
[forks-url]: https://github.com/nollive/projet-pacri/network/members
[stars-shield]: https://img.shields.io/github/stars/nollive/projet-pacri.svg?style=for-the-badge
[stars-url]: https://github.com/nollive/projet-pacri/stargazers
[license-shield]: https://img.shields.io/github/license/nollive/projet-pacri.svg?style=for-the-badge
[license-url]: https://github.com/nollive/projet-pacri/blob/main/COPYING
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://fr.linkedin.com/in/olivier-gaufr%C3%A8s