# Computational Neuroscience and Learning Models

MATLAB implementations of foundational models and algorithms in computational neuroscience and machine learning, completed as coursework at Sharif University of Technology.

## Topics Covered

**Neuron models**
- Hodgkin–Huxley model
- Izhikevich spiking neuron model
- Network of leaky integrate-and-fire neurons
- Noisy-output neuron model
- Bifurcation analysis and phase-plane plotting

**Spike-train analysis**
- Peri-Event Time Histogram (PETH) and raster plots
- Spike-Triggered Average (STA)

**Learning and networks**
- Hopfield network (auto-associative memory)
- Linear perceptron
- Two-layer neural network for handwritten-digit classification (MNIST-like)
- Reinforcement learning — classical conditioning paradigms
- K-Means clustering (image segmentation, includes a demo on `ponyo.jpeg`)
- Principal Component Analysis (1D toy data and an eigenfaces demo on `faces.mat`)

## Repository Structure

```
Bifurcation.m
Hodgkin-Huxley Model.m
Izhikevich Model.m
Linear Perceptron.m
Network of Integrate-and-fire Neurons.m
Noisy output Model.m
Phase Plain.m
Reinforcement Learning - Conditioning Paradigms.m
Hopfield Network/
K-Means/
PCA/
PETH and Raster Plot/
Spike Triggered Average/
two-layer neural network for hand digit problem/
```

## Requirements

- MATLAB (R2019b or later recommended)
- Statistics and Machine Learning Toolbox (for some scripts)
- Image Processing Toolbox (for the K-Means image demo)

## How to Run

Open MATLAB in this folder and run any script directly, e.g.:

```matlab
>> run('Hodgkin-Huxley Model.m')
>> cd('PCA'); run('PCA.m')
```

Subfolders contain their own `.mat` data files where needed.
