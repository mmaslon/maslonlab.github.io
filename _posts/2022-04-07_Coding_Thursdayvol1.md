
---
title: "Bash profile"
categories: [Coding Thursday]
---

We will be using shell to do most of our work. Shell accepts and runs commands. When you login, few files are read in to set up your environment. 

1) .bash_profile and .bashrc

.bash_profile is read first and it calls the second file, .bashrc. This is done every time you log in and it initializes and customizes your environment. 
.bash_profile file is where you can customize your environmental variable 
below a typical bash profile which also adds a directory ‘bin’ to your path ($PATH described above). 

<img src="https://github.com/mmaslon/maslonlab.github.io/blob/7cefb92e6420f29519b01073f0732e75e79db787/assets/theme/images/fig1_31333.png" width="400" align="left">
</figure>

2) $PATH 

a variable $PATH is a set of directories (separated by “ :” ) where the system assumes the executable programs are located. 

when you use an executable command in your script, the system will look into your $PATH to identify the executable code. 

3) bin

it is a good idea to create a sub-directory bin in your home directory - this is a place to keep any scripts or programs you write, 
or to add symlinks to some programs.  
