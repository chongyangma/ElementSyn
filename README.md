ElementSyn
==========

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/2854180796ac447a8eeeadea963dcda4)](https://app.codacy.com/gh/chongyangma/ElementSyn?utm_source=github.com&utm_medium=referral&utm_content=chongyangma/ElementSyn&utm_campaign=Badge_Grade_Settings)
[![Build Status](https://travis-ci.com/chongyangma/ElementSyn.svg?branch=master)](https://travis-ci.com/chongyangma/ElementSyn)
[![Build status](https://ci.appveyor.com/api/projects/status/p653ml95kprnmnn1?svg=true)](https://ci.appveyor.com/project/chongyangma/elementsyn)

This repository contains the source code and example data of the following publication:

> Dynamic Element Textures
>
> [Chongyang Ma](http://chongyangma.com/), [Li-Yi Wei](http://www.liyiwei.org/), [Sylvain Lefebvre](http://www.antexel.com/sylefeb/research), [Xin Tong](http://research.microsoft.com/en-us/um/people/xtong/xtong.html)
>
> In _ACM Transactions on Graphics (Proceedings of SIGGRAPH 2013)_
>
> [Project page](http://chongyangma.com/publications/dt/index.html),
> [Paper](http://chongyangma.com/publications/dt/2013_dt_paper.pdf),
> [Slides](http://chongyangma.com/publications/dt/2013_dt_slides.pdf),
> [Video](http://chongyangma.com/publications/dt/2013_dt_video.mp4),
> [YouTube](https://www.youtube.com/watch?v=dSvqGcBAorI),
> [BibTex](http://chongyangma.com/publications/dt/2013_dt_bib.txt)

Requirements
------------

The code requires [CMake](https://cmake.org/) to build and has been tested under Windows (MSVC 2010 and later), Linux and Mac OS X. Additional dependencies (included in this repo as submodules) are:
*   [eigen](https://gitlab.com/libeigen/eigen) for solving linear systems
*   [lodepng](https://github.com/lvandeve/lodepng) for saving sequences of screenshots

Usage
-----

Launching the compiled application from command line without any argument will print the usage information:

```bash
ParticleSystem.exe config_file.txt

TreeBranches.exe config_file.txt
```

Contact information
-------------------

Questions? Bug reports? Please send email to Chongyang Ma chongyangm@gmail.com .
