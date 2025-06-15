## CalVSP

CalVSP is a free, cross-platform, scientific software for volume, surface and PSA computations of single molecules and trajectory files. It is a command line oriented program. It can be executed in a shell (UNIX/Linux systems) or at a Windows or DOS command prompt.
![屏幕截图 2025-05-21 192116](https://github.com/user-attachments/assets/0ac588c7-f8f2-4dff-ae1c-41955d163c77)

## Software Requirement
If molecular conformation search, energy minimization of molecular volume, molecular surface area, and PSA calculation are required, it is recommended to install Open Babel and MOPAC.

Open Babel 3.0.0 -- Oct  7 2019 -- 20:03:12

MOPAC2016™


## **Installation**
The CalVSP program is written in C and has been tested on Linux and Windows platforms. It is provided as executable binary or source codes for each supported platform. To install the command line version of CalVSP please follow the instructions below.
1. Install GNU C compiler if you haven't already. On Ubuntu you may use sudo apt-get install gcc, for other Linux users you may need to use a different method.
2. Downloading the source code
3. compiling source code You can compile the program on your Linux computer using: 

```
gcc CalVSP.c -o CalVSP -lm -O3
```
  or

```
gcc CalVSP.c -o CalVSP -lm
```
If you are a Windows user, you can compile using gcc on the command prompt. Alternatively, install the free, portable, fast and simple C/C++ IDE Dev-C++, then Using Dev-C++, open the source code, click the compile option to compile and obtain the binary executable file.


## **How to Use**
Add CalVSP as an environment variable.
Open Terminal (Mac & Linux) or PowerShell (Windows).
The general synopsis for using CalVSP is:
```
CalVSP -i inputFile 
```
File extension of inputFile is very important (but case-sensitive) because CalVSP will judge an input file format according to its file extension. The MOL2, PDB, XYZ and SDF formats are supported. The XYZ format is recommended for trajectory files.

```
CalVSP -i alprenolol.xyz
```
If optional parameter -d is provided, CalVSP will also output surface points data. 

```
CalVSP -i inputFile -d suf_Data.XYZ
```

If optional parameter -s is provided, CalVSP will calculate solvent-accessible surface volume and PSA of molecules. If the optional parameter - s is not followed by radius data, the default water molecule radius of 1.4 Å will be used to calculate the solvent accessible surface area. Users can use the corresponding solvent molecule radius data to calculate the accessible surface area of the solvent based on the calculation situation.

```
CalVSP -i inputFile -s or CalVSP -i inputFile -s 1.4 or CalVSP -i inputFile -s 1.5
```

If the optional parameter - g is provided, CalVSP will perform a molecular conformation search on the input molecular structure data by calling the Obabel program. Then calculate the molecular surface area, molecular volume, and PSA data.

```
CalVSP -i inputFile -g
```

If the optional parameter - m is provided, CalVSP will optimize the molecular structure of the input molecular structure data by calling the Obabel program or MOPAC program. Then calculate the molecular surface area, molecular volume, and PSA data.
```
CalVSP -i inputFile -m
```
The selection parameters - g and - m can be used simultaneously to first search for molecular conformation and then optimize the molecular structure based on the input molecular structure data. Finally, the molecular surface area, molecular volume, and PSA data are calculated based on the optimized molecular data.
```
CalVSP -i inputFile -g -m
```



## **How to Cite**
## **License**
This package is distributed under the MIT License.It is provided with the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

<!--
**CalVSP/CalVSP** is a ✨ _special_ ✨ repository because its `README.md` (this file) appears on your GitHub profile.

Here are some ideas to get you started:

- 🔭 I’m currently working on ...
- 🌱 I’m currently learning ...
- 👯 I’m looking to collaborate on ...
- 🤔 I’m looking for help with ...
- 💬 Ask me about ...
- 📫 How to reach me: ...
- 😄 Pronouns: ...
- ⚡ Fun fact: ...
-->
