---
title: Interacting with the shell
---

## Interacting with the shell
<br>

A core function of the shell is to be the first running process a human user of
a computer can use to **interact with the kernel** to get things done. The Shell
achieves this by:

- accepting human oriented commands from a user.
- executing kernel system calls to get the work of the commands done.
- sending a response back. 

The shell includes the ability to locate other executable programs and launch new processes from them.
The shell is text based and is designed with programmers in mind. 
<br>

#### <ins>Shell variants</ins>

- Bourne
- BASH (Bourne Again SHell)
- C shell (csh)
- tcsh
- Korn Shell (ksh)
- Z Shell (zsh)

<br>


#### Command line

<br>
UNIX command lines can get quite complex. 

A hallmark of UNIX expertise is the ability to compose long command lines that chain together many commands. 

Learning to issue simple command lines is the first step to ruling the UNIX world.

When the shell receives a command line, it goes through a series of steps to process it. The rules of this processing define what is called the **Shell Grammar**.

Commands can be run ***interactively*** vs ***non-interactively (scripts)***


<br>
<left><img src="/img/commandline-flow.png" alt="commandline" width="100%"/></left>
<br>

#### Help commands

man pages  (manual pages)

Linux includes a built-in manual for nearly all commands, so these should be your go-to reference.

`man` 			format and display manual pages

1. NAME – a one-line description of what it does.
2. SYNOPSIS – basic syntax for the command line.
3. DESCRIPTION – describes the program's functionalities.
4. OPTIONS – lists command line options that are available for this program.
5. EXAMPLES – examples of some of the options available
6. SEE ALSO – list of related commands.

#### Shell command syntax
<br>
<left><img src="/img/shell-command-syntax.svg" alt="Shell command syntax" width="75%"/></left>
<br>