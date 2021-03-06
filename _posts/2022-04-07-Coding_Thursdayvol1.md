---
title: Bash profile
date: 2022-04-07 
categories: ["Coding Thursday"]
background: /assets/theme/images/bash.png
---

We will be using shell to do most of our work. Shell accepts and runs commands. When you login, few files are read in to set up your environment. 

1) .bash_profile and .bashrc

`.bash_profile` is read first and it calls the second file, `.bashrc`. This is done every time you log in and it initializes and customizes your environment. 
`.bash_profile` file is where you can customize your environmental variable 
below a typical bash profile which also adds a directory `bin` to your path (`$PATH` described above). 

```bash
# ~/.bash_profile: executed by bash(1) for login shells.
# see /usr/share/doc/bash/examples/startup-files for examples.
# the files are located in the bash-doc package.

# the default umask is set in /etc/login.defs
#umask 022

# include .bashrc if it exists
if [ -f ~/.bashrc ]; then
    . ~/.bashrc
fi

# set PATH so it includes user's private bin if it exists
if [ -d ~/bin ] ; then
    PATH=~/bin:"${PATH}"
fi
```

2) `$PATH` 

a variable `$PATH` is a set of directories (separated by `:` ) where the system assumes the executable programs are located. 

when you use an executable command in your script, the system will look into your `$PATH` to identify the executable code. 

to check what your path is you can use the following command:

```bash
echo $PATH
```

3) bin

it is a good idea to create a sub-directory `bin` in your home directory - this is a place to keep any scripts or programs you write, 
or to add symlinks to some programs.  

You should add bin to your bash profile [it is already done in the example above].


