# FBG-Simulation-PyQt

FBG Simulation PyQt is an open-source project that builds upon the previous work
of [Gilmar Pereira](https://github.com/GilmarPereira)'s
project [FBG_SiMul](https://github.com/GilmarPereira/FBG_SiMul) and
[Benjamin Frey](https://github.com/benfrey)'s project
[FBG-SimPlus](https://github.com/benfrey/FBG-SimPlus).

Our goal is to provide a clean sample implementation of a Qt application that
can be used in any simulations projects. It proves how to effectively decouple
simulation code from GUI which is really tough to manage without years of
experience in developing Qt or GUI.

And hopefully allows scientist and researchers to write richer simulation
without the hassle of writing a stable and clean GUI application.

This project is released under [BSD 3](./LICENSE), and we adhere to the
licensing and attribution requirements of the original projects that we have
adapted.

## Features

- Easy-to-use interface for managing simulations
- Decoupled implementation of FBG simulation (fully tested with pytest)
- [WIP, will add more]


## Getting Started

### Prerequisites

The python application is built with [PySide6](https://pypi.org/project/PySide6/)
and uses [pytest](https://docs.pytest.org/) for all its unit tests.

You can install all required packages with the followin command:

```bash
pip install -r requirements.txt
```


### Installation

A big WIP!

### Usage

To start the application simply run the [main.py](./main.py) and load the
provided [sample.txt](./sample/tut-export.txt) datafile.


## Acknowledgments and Licensing

This project is built upon, and adapted from, the following open source projects:

- [FBG_SiMul](https://github.com/GilmarPereira/FBG_SiMul), licensed under
  [GPL-3](https://github.com/GilmarPereira/FBG_SiMul/blob/master/LICENSE)
- [FBG-SimPlus](https://github.com/benfrey/FBG-SimPlus), licensed under
  [GPL-3](https://github.com/benfrey/FBG-SimPlus/blob/master/LICENSE)

Please make sure to review the licenses of these projects to understand their terms and conditions. Our project adheres to these licenses and provides proper attribution to the original authors.

**Citation**: If you choose to use elements of this application in your own
work, please cite the following paper:

> Frey, B., Snyder, P., Ziock, K., & Passian, A. (2021).
> _"Semicomputational calculation of Bragg shift in stratified materials"_.
> Physical Review E, 104(5), 055307.


## License

This project is licensed under the BSD-3 license. Please read the
[LICENSE](./LICENSE) file for full details.
