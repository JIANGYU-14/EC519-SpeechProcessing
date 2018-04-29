# Speech Processing by Humans and Machines Project

## Porject Introduction

The project is based on the non-quantized system, forwardly, vocal track coefficients and pitch value has been quantized to exhibit better synthesizing performance on given speech file, binary impulse generator included in the synthesizing process been used for produce excitation, and residual error represented via pitch pulses and noise(binary) through LPC coding.

In the second part, CELP(Code Excited Linear Prediction) algorithmhas been implemented for speech file synthesizing, similarly, the system should be quantized as well. The algorithm is based on fourmain ideas as below.

* Using the source-filter model of speech production through linear prediction (LP)
* Using an adaptive and a fixed codebook as the input (excitation) of the LP model
* Performing a search in closed-loop in a â˘AIJperceptually weighted domainâ
* Applying vector quantization (VQ)

In a nutshell, the project is organized in a way that binary source which served as excitation should first be tested, and effect of quantization on synthesized speech quality illustrated in figures. Followed by that, CELP coding has been used for speech synthesizing, and the residual represented by codewords from a VQ-generated codebook after long-term(pitch period) and short-term(vocal tract) prediction on each frame, rather than by multiple pulses. Finally, all the parameters follow the same value as in previous project.

## Binary Source

### Block Diagram

### Parameter Set
Parameters involved in the project are displayed as below
* HammingWindow Length =Window Shift = 30ms (240Samples)(Nonoverlapping frames)
* Parameter Set: Log Area Ratio Coefficients (As spectral sensitivity measures sensitivity to errors in the log area ratio parameters, considering the bit-rate budget, this set of parameters much less sensitive than others)
* Prediction Order P = 10 Specified in Requirement
* Bit allocation: 7 7 7 6 6 5 5 4 4 3; Gain allocation: 6 (2400 bps)
* Bit allocation: 12 10 10 9 9 8 8 8 7 7; Gain allocation: 8 (3600 bps)
* Bit allocation: 16 16 15 14 14 14 13 13 11 9; Gain allocation: 13 (4800 bps)
* Bit allocated to each position decreasing as order increasing, which comes along with the decreasing of short term energy, and bit allocated to gain is approximately at the mean of LPC-10 bit allocation

## CELP(Code Excited Linear Prediction)

### Block Diagram

### Parameter Set

Parameters follow the same value as before, and implemented on 4800bps and 9600 bps.



