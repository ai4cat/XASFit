# **XANESFit-Screen: Automated XANES Peak Fitting and FEFF10-Based Configuration Screening**

An automated Python workflow for XANES spectral analysis, peak fitting, and structure screening by combining experimental spectra with FEFF10-simulated references.

This repository provides a reproducible pipeline for:

- automated peak detection in experimental XANES spectra,
- Gaussian/Voigt-based spectral fitting,
- descriptor extraction from fitted spectra,
- comparison between experimental and FEFF10-simulated spectra,
- ranking and screening of candidate local configurations.

The workflow is designed to accelerate XANES interpretation while reducing manual intervention in peak assignment and spectrum matching.

For further details, you are encouraged to consult our [paper](https://) and visit our [website](https://) for additional resources and also our [dataset](http://openadc.com.cn:23345/) information.

## **Installation**

### Development Environment
- Python 3.13
- Validated on Linux OS
- Use `conda env create -f xasfit.yml` to create the enviornment.

### Setup
To set up the codes, run the following commands:

```bash
git clone https://github.com/ai4cat/XASFit.git
cd XASFit
```

## Data Availability & Copyright Notice

The experimental XANES spectra used in this project were obtained through in-house measurements conducted within our research group. The corresponding atomic configurations and theoretical models were independently constructed as part of our internal database development.

Due to data ownership and ongoing research considerations, the full dataset (including experimental spectra and computed structural configurations) is not publicly distributed within this repository.

Interested users may access relevant data and explore the curated database through our official platform:

👉 [Open-ADC](http://openadc.com.cn:23345/)

For specific requests or collaborations, please feel free to contact us.

## **Repository Structure**

```text
XASFit/
├── fitting/
└── pre-process/                        
```

## Contributing
Contributions are welcome! Please follow the standard fork-and-pull request workflow on GitHub.

If you use our code in your research, please cite our paper:
```bash
@article{,
  title={s},
  author={},
  journal={},
  year={},
  volume = {},
  pages = {}
}
```

## **License**

This project is licensed under the [CC-BY-ND-NC License](https://github.com/ai4cat/XASFit/blob/main/LICENSE). Please see the LICENSE file for more details.

