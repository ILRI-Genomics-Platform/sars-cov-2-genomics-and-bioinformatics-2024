---
title: Filesystems in UNIX
---

## Filesystem in UNIX/LINUX
<br>

A ***filesystem*** can be defined as a data structure or a collection of files. 

UNIX system regards everything as a file. Files can be divided into three categories; 

- Ordinary or plain files.
- Directories.
- Special or device files.

In UNIX, filesystem can refer to two very distinct things; the **directory tree** or the **arrangement of files on disk partitions**. 
The latter can be thought of as the physical filesystem as it has a tangible physical location.

<br>


#### Filesystem: files

A file is an abstract object that can store information. We can write and read contents to and from it. 

***Data*** of a file contain a collection of bytes. The bytes of the file might encode ASCII data, image data, audio data, an executable program, etc.

***Metadata*** of a file are associated descriptive facts. Examples of these facts include:
- who owns the file
- the length of the file (measured in bytes)
- who has permissions to read or write its contents
- the time the contents was last modified
- the time that it was last read
- the time that the descriptive facts where last changed (e.g., the file
  permissions were modified)
<br>
<left><img src="/img/filesystem-files.png" alt="filesystem-file" width="75%"/></left>
<br>


#### Filesystem: directory tree

<br>
<left><img src="/img/filesystem-directory-01.png" alt="filesystem-directory" width="100%"/></left>
<br>

<br>
<left><img src="/img/filesystem-directory-02.jpg" alt="filesystem-directory" width="75%"/></left>
<br>

#### Paths
In UNIX, the top of the tree is the one directory that always exists and is not a sub-directory of any other directory. This is the `root` directory.


The name of the variables file at the bottom of the diagram is composed by joining the name of the directories leading to it along with its name at the end.

`/ + home + jolo + Project`

 `/ + home + jolo + Results`

  `/ + home + jdoe + Documents`




Notation for a full path name of a file or directory is by joining the
independent components with the “/”.
```bash
    /home/jolo/Project
    /home/jolo/Results
    /home/jdoe/Documents
```


Root path `/`

<br>


## <ins>Absolute paths</ins>

Details the entire path through the directory structure to get to a file, starting at /, the root directory

```/home/jolo/Projects```

## <ins>Relative paths </ins>

Is the path from where you are now (your present working directory) to the file/directory in question

```../Results```


<br>

## QUIZ: Test you knowledge on UNIX/LINUX filesystem structure
<br>

<br>
<left><img src="/img/filesystem-directory-03.png" alt="filesystem-directory" width="100%"/></left>
<br>

1. A UNIX file is a physical object that exists on the computer’s hard-drive: True or False?
2. What is meta-data versus data with respect to files?
3. The directory structure and name of the directories is fixed: True or False?
4. In the directory tree example in the lecture what is the full path name of user01’s workflows directory?
5. In the directory tree example in the lecture what is the full path name of user01’s assignment files?
A directory can contain both files and directories: True or False?
<br>
