# Configuring PDAF in EAT

The straight forward solution would have been including the necessary  CMakeLists.txt files in the PDAF Github repository - but this has been  refused by PDAF/Lars Nerger with the explanation that it would risk jeopardising the default Make based build system. Despite a lengthy email correspondance with several example configuration command lines no agreement could be made.

To circumvent this obstacle it has been necessary to include support for building PDAF as an integrated part of building EAT without having the necessary CMakeLists.txt files inside the PDAF directory structure. Unfortunately this method is a bit more fragile than if the CMake configuration was an integrated  part of PDAF. The issue being that there is not an easy way to guarantee that the version of PDAF pointed to by the PDAF_BASE configuration variable and the version of PDAF specified during EAT configuration corresponds.

Anyway here is the recipe:

1) Not specifying any -DPDAF_BASE during configuration will use the 
git submodule PDAF version checked out as part of a standard EAT cloning process. No further considerations are necessary.

2) When specifying -DPDAF_BASE it must be ensured that the PDAF configuration version pointed to by PDAF_BASE is corresponds to what is in .../pdaf/CMakeLists.txt.

The second point is mainly for developers and unless there is a good reason just configure with default settings.

It is possible to compile an executable that will show some configuration details. It works like this (using make as an example).
In the build directory:

```
make -j 19 pdaf_configure
```

If no errors:
```
pdaf/v2.1/src/pdaf_configure 
```

![image-20230427134500330](/home/kb/.config/Typora/typora-user-images/image-20230427134500330.png)
