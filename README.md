## CalVSP

CalVSP is a free, cross-platform, scientific software for volume, surface and PSA computations of single molecules and trajectory files. It is a command line oriented program. It can be executed in a shell (UNIX/Linux systems) or at a Windows or DOS command prompt.
![屏幕截图 2025-05-21 192116](https://github.com/user-attachments/assets/0ac588c7-f8f2-4dff-ae1c-41955d163c77)


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
```
CalVSP -i alprenolol.xyz
```

## **How to Cite**
## **License**
This package is distributed under the MIT License.

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
