# MWF-MIMOSA

MWF-MIMOSA is a Pulseq-based MRI acquisition and reconstruction framework for myelin water fraction mapping.

Please run `setup.m` first to initialize the required MATLAB paths.

## Pulseq Sequence

- Download the Pulseq sequence programming environment from [Pulseq](https://pulseq.github.io/) and add it to the MATLAB path.
- Run `01_gen_Pulse_Seq/write_MWF_MIMOSA_1mm.m` to generate the Pulseq sequence file for the 1 mm isotropic protocol.
- Run `01_gen_Pulse_Seq/write_MWF_MIMOSA_700um.m` to generate the Pulseq sequence file for the 700 µm isotropic protocol.

## Reconstruction

The baseline zero-shot reconstruction code is forked from:

- [ZS-SSL](https://github.com/byaman14/ZS-SSL)
- [Zero-DeepSub](https://github.com/yohan-jun/Zero-DeepSub)

### Installation

The required dependencies are provided in:

```bash
02_Recon/zsssl_recon_1mm/environment_tf2.yml
```

They can be installed using:

```bash
conda env create -f 02_Recon/zsssl_recon_1mm/environment_tf2.yml
```

### Data

The raw MWF-MIMOSA data with 1 mm isotropic resolution can be downloaded from [here](https://drive.google.com/file/d/1KBWsIbsEYRE9UPjvUFiUG7LNIos81G2t/view?usp=drive_link).

After downloading the raw data, place it in:

```bash
02_Recon/rawdata
```

### Reconstruction Pipeline

#### 1. Preprocessing

Run the following MATLAB script to prepare the raw data for zero-shot self-supervised reconstruction:

```matlab
02_Recon/prepare_data_for_zsssl_recon.m
```

#### 2. Training

Run one of the following scripts to perform multi-contrast and multi-slice zero-shot self-supervised learning training:

```bash
python 02_Recon/zsssl_recon_1mm/zs_ssl_train_multi_mask_batch_v10_ms.py
```

or

```bash
python 02_Recon/zsssl_recon_700um/zs_ssl_train_multi_mask_batch_v10_ms.py
```

Before training, reconstruction hyperparameters can be adjusted in `parser_ops.py` under the corresponding reconstruction folder.

#### 3. Inference

Run one of the following notebooks to load the checkpoints saved during training and perform inference:

```bash
02_Recon/zsssl_recon_1mm/zs_ssl_inference_ms.ipynb
```

or

```bash
02_Recon/zsssl_recon_700um/zs_ssl_inference_ms.ipynb
```

## Parameter Estimation

The parameter estimation pipeline is adapted from:

- [GACELLE](https://github.com/kschan0214/gacelle)
- [MIMOSA](https://github.com/yutingchen11/MIMOSA)
- [chi-separation](https://github.com/SNU-LIST/chi-separation)

### 1. Dictionary Generation

Run one of the following MATLAB scripts to generate the dictionary:

```matlab
03_ParamEstimation/gen_MWF_MIMOSA_dict_1mm.m
```

or

```matlab
03_ParamEstimation/gen_MWF_MIMOSA_dict_700um.m
```

### 2. Mapping

Run one of the following MATLAB scripts to perform parameter estimation:

```matlab
03_ParamEstimation/mapping_mwf_mimosa_1iso.m
```

or

```matlab
03_ParamEstimation/mapping_mwf_mimosa_700um.m
```

## Copyright and License Notice

This project is licensed for non-commercial research use only.

For other purposes, please contact:

```text
ychen156@mgh.harvard.edu
```
