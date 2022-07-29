# DirectTranscription

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://GrantHecht.github.io/DirectTranscription.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://GrantHecht.github.io/DirectTranscription.jl/dev)
[![Build Status](https://github.com/GrantHecht/DirectTranscription.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/GrantHecht/DirectTranscription.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/GrantHecht/DirectTranscription.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/GrantHecht/DirectTranscription.jl)

A Julia package for solving multi-phase optimal control problems with direct transcription. DirectTranscription.jl currently supports several Lobatto IIIA collocation algorithms and Radau pseudospectral collocation is in the works. The resultant NLP problem can then be solved using Ipopt (default) or Snopt (requires licence) through the Snopt.jl package developed by the BYU Flow Lab (see [link](https://github.com/byuflowlab/Snopt.jl)). To use Snopt, the Snopt.jl source code must be cloned and built as described at the previously provided link.

Currently, no documentation is written for DirectTranscription.jl but this is in the works as well. In the meantime, the examples folder contains several scripts for solving the Brachistichrone problem and a low-thrust transfer between halo orbits in the Earth-Moon system. Feel free to email me with any questions at granthec@buffalo.edu.
