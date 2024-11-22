# **Meraxes**
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit) [![All Contributors](https://img.shields.io/badge/all_contributors-2-orange.svg?style=flat-square)](https://github.com/meraxes-devs/meraxes/graphs/contributors)  

<figure>
  <img src="https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/mnras/526/1/10.1093_mnras_stad2448/1/m_stad2448fig2.jpeg?Expires=1734589748&Signature=P0bapPBEJaZqeETwWL7Uo0ES7VDOEov1--cFCO2oGhwVQRCWz8RRnDA9~42DViUVMoAQJChkkFX5AwgxY6oQVq7AY8MZjrls4LAdOR-EjVZixnhufCh5-fo11-M4MOONglvOZNeuu2GydKswOzgc6qqGyJb0u3eHNjnfVcAiefDZtR5bPp6ZFDycBtW9nYr05HPL5NSUP4RUnMFqOLX5aNDIPyNBktNGuxbmrf~Id1HOsIz07NgQqBHRZDZ383GUe-kFWX1ydFjcpF1X--7ZBdjZh-zUYamRS6pJGCxDwOXT7dOhlV8uxiRcSKd4W2enzXzEFu-gG6yRUBtztF3TVw__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA" style="width: 100%; height: auto;" alt="Galaxy UV luminosity function"/>
  <figcaption>Galaxy UV non-ionizing luminosity functions from z = 20 to 5 predicted by our fiducial model.</figcaption>
</figure>


<figure>
  <img src="https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/mnras/520/3/10.1093_mnras_stad281/1/m_stad281fig7.jpeg?Expires=1735221664&Signature=eCtJmKVKts0UGgNKMm75Sx6VqgFkY14abE7ifwWCCBTa~VAgvstvcH03-AiMgalY4GzQXJVY8tckN47TZTJY9iFMZnEoT1oZKg3rgUKQUSJoYhEs1U8p8BLRk9MypDzUegeYBvMCIU7dnvGFE2NbwABahHCJWrhiYEAlpQ-ICd-zP~h7dgMs1DGnX3N~NB-bUlZd5A4RM1Q7Ep3LjRydM719MQG3mU0OsirPHJSDvAtsTpg-9Pg-4VW-QnnapyA2slCmYklv9KHwkIgBcLJTCVbL2gpuOBGD7FdfZa1~9xbdqK~9Xakhv~fkAL-QG-IyA--RMJbs8vVFuya0bS71aw__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA" style="width: 100%; height: auto;" alt="Galaxy UV luminosity function"/>
  <figcaption>The light-cone evolution of the 21-cm brightness temperature from our simulations.</figcaption>
</figure>



Meraxes is a **semi-analytic galaxy formation model** designed to study the interplay between high-redshift galaxies and their intergalactic medium (IGM) during the Epoch of Reionization (EoR). It  evaluates galaxy properties (e.g., stellar mass, star formation rate, UV luminosity in rest frame or through observational filters such as HST, JWST or other user-defined ones) from dark matter halo merger trees. It computes the photon budget including UV and X-ray, integrating with a customized version of [``21cmFAST``](https://github.com/21cmfast/21cmFAST) to simulate the ionization state and 21-cm signals of the IGM during the EoR.  

A full list of publications directly using Meraxes is available on [ADS](https://ui.adsabs.harvard.edu/public-libraries/CWUcYnt3TsmG6BuOKjR0Fw).  

## **Complementary Tools**
- **[DRAGONS](https://github.com/meraxes-devs/dragons):** A Python package for reading and processing Meraxes output.  
- **[Sector](https://github.com/meraxes-devs/sector):** A library for computing spectral energy distributions from Meraxes simulations.  


## **Installation**
Please refer to [`BUILD.md`](./BUILD.md) for detailed installation and build instructions.


## **Documentation**
Comprehensive documentation is currently under development. For inquiries or assistance, please contact the team.  


## **Acknowledging**
If you use Meraxes or its outputs in your research, please cite the following foundational paper:  

- Mutch et al. (2016). **Dark-ages reionization and galaxy formation simulation - III. Modelling galaxy formation and the epoch of reionization.**  *Monthly Notices of the Royal Astronomical Society*, 462 (1): 250–276. [DOI: 10.1093/mnras/stw1506](https://doi.org/10.1093/mnras/stw1506).  

If using specific features introduced in Meraxes, please cite the corresponding papers:  

- **AGN Model:**  Qin et al. (2017). **Dark-ages reionization and galaxy formation simulation - X. The small contribution of quasars to reionization.** *Monthly Notices of the Royal Astronomical Society*, 472 (2): 2009–2027. [DOI: 10.1093/mnras/stx1909](https://doi.org/10.1093/mnras/stx1909).  

- **SED and Dust Model:**  Qiu et al. (2019). **Dark-age reionization and galaxy formation simulation - XIX. Predictions of infrared excess and cosmic star formation rate density from UV observations** *Monthly Notices of the Royal Astronomical Society*, 489 (1): 1357-1372. [DOI: 10.1093/mnras/stz2233](https://doi.org/10.1093/mnras/stz2233).  

- **21-cm Model:**  Balu et al. (2023). **Thermal and reionization history within a large-volume semi-analytic galaxy formation simulation.**  *Monthly Notices of the Royal Astronomical Society*, 520 (3): 3368–3382. [DOI: 10.1093/mnras/stad281](https://doi.org/10.1093/mnras/stad281).  

- **Minihalo Model:**  Ventura et al. (2024). **Semi-analytic modelling of Pop. III star formation and metallicity evolution - I. Impact on the UV luminosity functions at z = 9–16.** *Monthly Notices of the Royal Astronomical Society*, 529 (1): 628–646. [DOI: 10.1093/mnras/stae567](https://doi.org/10.1093/mnras/stae567).  

## **Contributors**

<table>
  <tbody>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://smutch.github.io/"><img src="https://avatars.githubusercontent.com/u/782987?v=4?s=100" width="100px;" alt="Simon Mutch"/><br /><sub><b>Simon Mutch</b></sub></a><br /><a href="https://github.com/meraxes-devs/meraxes/commits?author=smutch" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://researchportalplus.anu.edu.au/en/persons/yuxiang-qin"><img src="https://avatars.githubusercontent.com/u/15994713?v=4?s=100" width="100px;" alt="Yuxiang Qin"/><br /><sub><b>Yuxiang Qin</b></sub></a><br /><a href="https://github.com/meraxes-devs/meraxes/commits?author=qyx268" title="Code">💻</a> <a href="https://github.com/meraxes-devs/meraxes/pulls?q=is%3Apr+reviewed-by%3Aqyx268" title="Reviewed Pull Requests">👀</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/EMventura"><img src="https://avatars.githubusercontent.com/u/98299102?v=4?s=100" width="100px;" alt="Emanuele Maria Ventura"/><br /><sub><b>Emanuele Maria Ventura</b></sub></a><br /><a href="https://github.com/meraxes-devs/meraxes/commits?author=EMventura" title="Code">💻</a> <a href="https://github.com/meraxes-devs/meraxes/pulls?q=is%3Apr+reviewed-by%3AEMventura" title="Reviewed Pull Requests">👀</a></td>
      <td align="center" valign="top" width="14.28%"><a href="http://s-balu.github.io"><img src="https://avatars.githubusercontent.com/u/14290533?v=4?s=100" width="100px;" alt="Balu Sreedhar"/><br /><sub><b>Balu Sreedhar</b></sub></a><br /><a href="https://github.com/meraxes-devs/meraxes/commits?author=s-balu" title="Code">💻</a> <a href="https://github.com/meraxes-devs/meraxes/pulls?q=is%3Apr+reviewed-by%3As-balu" title="Reviewed Pull Requests">👀</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/BradGreig"><img src="https://avatars.githubusercontent.com/u/16087482?v=4?s=100" width="100px;" alt="BradGreig"/><br /><sub><b>BradGreig</b></sub></a><br /><a href="https://github.com/meraxes-devs/meraxes/commits?author=BradGreig" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="http://www.astronomy.swin.edu.au/~gpoole/"><img src="https://avatars.githubusercontent.com/u/599836?v=4?s=100" width="100px;" alt="Gregory B. Poole"/><br /><sub><b>Gregory B. Poole</b></sub></a><br /><a href="https://github.com/meraxes-devs/meraxes/commits?author=gbpoole" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/pgeil"><img src="https://avatars.githubusercontent.com/u/13758421?v=4?s=100" width="100px;" alt="pgeil"/><br /><sub><b>pgeil</b></sub></a><br /><a href="https://github.com/meraxes-devs/meraxes/commits?author=pgeil" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/yqiuu"><img src="https://avatars.githubusercontent.com/u/26683739?v=4?s=100" width="100px;" alt="Yisheng Qiu"/><br /><sub><b>Yisheng Qiu</b></sub></a><br /><a href="https://github.com/meraxes-devs/meraxes/commits?author=yqiuu" title="Code">💻</a></td>
    </tr>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/daviesje"><img src="https://avatars.githubusercontent.com/u/36873665?v=4?s=100" width="100px;" alt="daviesje"/><br /><sub><b>daviesje</b></sub></a><br /><a href="https://github.com/meraxes-devs/meraxes/commits?author=daviesje" title="Code">💻</a></td>
    </tr>
  </tbody>
</table>
