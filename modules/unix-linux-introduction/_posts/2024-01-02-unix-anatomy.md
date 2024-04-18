---
title: What is UNIX
---

### What is UNIX
<br>
UNIX is an operating system (OS) first developed in the 1960s by Ken Thompson and Dennis Ritchie. 
Still under constant development ever since. This operating systems are widely
used in scientific computing (including bioinformatics)
<br>

#### Features

1. Stable, multi-user, multi-tasking system for servers, desktops and laptops.
2. Flexible and easy to adopt to specific needs. UNIX philosophy – ‘small is good’ - each program is designed to do one job well.
3. Written in a machine independent language – UNIX can run on various hardware (laptops, desktops, super-computers).
4. Open source/freely available with open-source codes.
5. UNIX systems also have a graphical user interface (GUI) like Microsoft
   Windows which provides an easy-to-use environment. 
<br>

#### Unix anatomy
<br>

UNIX is multi-layered system

<left><img src="/img/unix-anatomy.png" alt="The anatomy of a UNIX system" width="75%"/></left>
The anatomy of a UNIX system.

<br>

#### <ins>Kernel</ins>

The kernel is the central program in an operating system and has ultimate authority over the computer.
Bottom layer of the system that provides a means for starting application programs (often also called user programs).

The kernel facilitates four basic types of services:
- creation and management of processes a filesystem
- communications
- a means to start the system

Kernel functions, such as **allocation of memory** and **CPU**, are performed without being explicitly requested by user processes.

Other functions of the kernel, such as **resource allocation** and **process creation and management**, are initiated by **requests from processes**.

These requests from processes come in the form of **system calls**. A system call can be thought of as a low-level request to the operating system (e.g., kill, fork, exec, exit, close)
<br>

#### <ins>Process</ins>

In UNIX, a running program is called a process. A program that a process can be created from is called and executable.

#### <ins>The Shell</ins>

UNIX system user accesses the services of the kernel through an interface called a shell. 

The shell is a command interpreter that allows the user to initiate processes to perform a nearly infinite variety of tasks.


#### <ins>Utilities</ins>

Utilities are programs that perform system functions. Utility can also refer to a command that is used to do work of some sort, such as mv to move files or directories.


## Linux

A Linux distribution, often known as a Linux distro, is a whole operating system
based on the Linux kernel. It is made up of the Linux kernel, core system libraries,
user applications, utilities, and a variety of software packages that work together to
offer a functional computer environment.
<br>


#### Why is LINUX commonly used by bioinformaticians

- Linux is commonly used in scientific computing because it is Open Source and has
attracted many different contributors.
- Linux's command line interfaces CLI enables bioinformaticians to handle data, automate
tasks, and run complex research pipelines. The CLI processes large volumes of data, text
files, and bioinformatics applications.
- Linux-optimized bioinformatics tools are plenty. Alignment, sequence analysis, genomics,
proteomics, and statistical analysis methods are used by researchers.
- Scalability and speed allows Linux to perform computationally difficult tasks including
Next Generation Sequencing (NGS) and genomic analysis. Bioinformaticians can analyze
data quicker using multi-core processors and distributed computing systems
<br>

#### Commonly used Linux Distributions
1. Ubuntu: Ubuntu is one of the most widely used Linux distributions. It focuses on ease of
use and provides a user-friendly interface. It has both a desktop version and a server
version.
2. Debian: Debian is a stable and highly customizable distribution. It emphasizes free and
open-source software. Many other Linux distributions, including Ubuntu, are derived from
Debian.
3. Fedora: Fedora is a community-driven distribution sponsored by Red Hat. It aims to
showcase the latest technologies and software in the Linux ecosystem. It is known for its
strong focus on security and new features.
<br>

