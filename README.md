<!-- PROJECT LOGO -->
<br />
<p align="center">
  
  <img src= https://github.com/martindoff/Energy-management-hybrid-aircraft/blob/main/Pic1.jpg alt="Logo" width="600" height="300">
  <p align="center">
   Predictive energy management for hybrid-electric aircraft with parallel propulsion system using convex optimisation.
    <br />  
  </p>
</p>



<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

The following case study considers the problem of computing the optimal power split for a hybrid-electric aircraft with a parallel propulsion system (gas turbine and electric motor) using convex optimisation. The goal is to obtain the optimal ratio between electric and gas turbine usage in order to meet the drive power demand. We consider solving the problem with 1) the generic convex programming package CVX, 2) the Alternating Direction Method of Multipliers (ADMM), 3) a heuristics (CDCS).

The theory behind the algorithms is presented in the paper: "Predictive energy management for hybrid electric aircraft propulsion systems" by Martin Doff-Sotta, Mark Cannon and Marko Bacic. 

### Built With

* MATLAB
* CVX
* Mosek



<!-- GETTING STARTED -->
## Getting Started


### Prerequisites

You need to install CVX for MATLAB:

* [cvx](http://cvxr.com/cvx/download/)

### Instructions 

1. Clone the repository
   ```sh
   git clone https://github.com/martindoff/Energy-management-hybrid-aircraft.git
   ```
2. Copy the directory into your MATLAB path
3. Run the script
   ```matlab
   main
   ```

<!-- ROADMAP -->
## Roadmap


<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<!-- CONTACT -->
## Contact

Martin Doff-Sotta - martin.doff-sotta@eng.ox.ac.uk

Linkedin: https://www.linkedin.com/in/mdoffsotta/



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/github_username/repo.svg?style=for-the-badge
[contributors-url]: https://github.com/github_username/repo/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/github_username/repo.svg?style=for-the-badge
[forks-url]: https://github.com/github_username/repo/network/members
[stars-shield]: https://img.shields.io/github/stars/github_username/repo.svg?style=for-the-badge
[stars-url]: https://github.com/github_username/repo/stargazers
[issues-shield]: https://img.shields.io/github/issues/github_username/repo.svg?style=for-the-badge
[issues-url]: https://github.com/github_username/repo/issues
[license-shield]: https://img.shields.io/github/license/github_username/repo.svg?style=for-the-badge
[license-url]: https://github.com/github_username/repo/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/github_username
