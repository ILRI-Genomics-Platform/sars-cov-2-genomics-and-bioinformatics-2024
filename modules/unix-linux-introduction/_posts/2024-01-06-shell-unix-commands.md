---
title: UNIX shell commands
---

## UNIX shell commands
<br>

#### <ins>Listing files and directories</ins>

`ls` 			lists the contents of your current working directory.
`ls -a` 	 	lists all the contents of your current working directory (including directories with a dot)

**Exercise**
1. What does the command `ls –lth` do?
2. To list all the files in **long** format and show complete time information, which command can you use?
3. To list all the files in a directory and output in a **comma-separated format**, which optional flag can you use?
<br>


#### <ins>Making directories</ins>

`mkdir` 			make a directory.

**Exercise**
Make a directory in your `home` directory and name it `lectures`
In `lectures` directory make additional sub-directories, `unix-stuff` and `shell-scripting`
Make another directory in your `home` directory and name it `assignments`



`cd`    change directory.
`cd .`  current directory
`cd ..` parent directory

Typing `cd` with no argument always returns you to your `home` directory.
<br>


#### <ins>Retrieve data/files over network</ins>

`wget`              non-interactively download of files from the Web.  
                    It supports HTTP, HTTPS, and FTP protocols as well as 
                    retrieval through HTTP proxies

**Exercise**
Download the `science.txt` file to your `home` directory.

`wget https://hpc.ilri.cgiar.org/~jjuma/training_data/science.txt`
<br>

#### <ins>Current working directory</ins>

`pwd`               write absolute pathname of the current working directory to the standard output


#### <ins>Copying files</ins>

`cp`	 			copy files

Typical command

`cp file1.txt file2.txt`		copy `file1.txt` in the current working directory and name it `file2.txt`

**Exercise**
Copy the `science.txt` file to the `assignments` directory.
Assuming you have a file named file1.txt in the current working directory, what
does the command `cp file1.txt .` do?
<br>


#### <ins>Moving files</ins>


`mv`	 			move or rename files

Typical command
`mv source target`		    rename 
`mv source directory`		move

**Exercise**
In your `home` directory, make a directory named `backups` and move the `science.txt` file from assignments to the `backups` directory.
<br>


#### <ins>Clearing the terminal screen</ins> 


`clear`	 			clear the terminal screen

Display content of a file

`cat`	 			concatenate and print files
`less`	 			Allows scrolling and searching through output.
`more`	 			display contents of a text file on terminal screen one page at a time
`head`	 			display first lines of a file
`tail`	 			display last lines of a file

**Exercise**
- Display all the contents of the `science.txt` file
- Display the **last 5 lines** in the `science.txt` file
- Display the **first 5 lines** in the `science.txt` file
<br>


#### <ins>Searching the contents of a file</ins>


`less`	 			display contents of a text file on terminal screen one page at a time
To search for a word in a text file, we can use the less command followed by the forward slash (/) and type the word to search. This will highlight all the occurrences of the word. Type [n] to search for the next occurrence of the word.

**Exercise**
Search for all the occurrences of the word science in the `science.txt` file
<br>


#### <ins>Searching the contents of a file with grep (global regular expression print)</ins>


`grep`	 			searching and matching patterns in a file
`grep` is case sensitive. You can ignore this feature by using the flag `–i`, i.e `grep –i`

**Exercise**

`grep science science.txt`

- Search for all occurrences of the word **science**, ignore case sensitivity.
- Display lines that **do not match the word science**.
- Print only the **total count** of the matched lines.
- Display the **line that match preceded by the line number**.
<br>


#### <ins>Word and line counts</ins>


`wc`	 			word, line, character and byte count


**Exercise**
`wc science.txt`

- How many **lines** are in the `science.txt` file.
- How many **characters** are in the science.txt file.
<br>

#### <ins>I/O and Redirection</ins>
<br>
Most processes in UNIX:

standard input (`stdin`): stream data going into a program. By default, this is input from the keyboard.
standard output  (`stdout`): output stream where data is written out by a program. By default, this output is sent to the screen.
standard error (`stderr`), another output stream (independent of stdout) where programs output error messages. By default, error output is sent to the screen.

The output of a process or program can be saved using the redirection operator `>`.
For example, you can redirect the output of the command ls to a file by:
    `ls –lth > list_of_files.txt`
This will create a file called `list_of_files.txt` or overwrite it if it does not exist.

If you know that the file exists, you can append more content by using the
operator `>>`
    `ls –lth >> list_of_files.txt`
<br>



Redirect the standard input into a file using > operator.

`cat > list1`

Type the following in an interactive manner
pear
banana
apple
Hit `Ctrl ^ [d]` to stop the input stream

Using the above method, create another file called `list2` containing the following fruits: orange, plum, mango, grapefruit. Read the contents of `list2`

Concatenate the two lists, `list1` and `list2` and name the resulting list as
`biglist`.

#### <ins>Sort items</ins>

`sort biglist`


#### <ins>Input Redirection</ins>

Input can also be given to a command from a file instead of typing it in the shell by using the redirection operator `<`.
`cat < science.txt`


#### <ins>Error Redirection</ins>

`rmdir no-dir 2> nodir_error.txt`

#### <ins>Merge stdout with stderr</ins>

`cat science > combined_output.txt 2>&1`

**Exercise**
Use `biglist` as input stream, `sort` the list and **redirect** to a standard
output list called **slist**

<details markdown="1">
<summary>Show answer</summary>

```bash
cat < biglist | sort > slist
```
</details>
<br>


#### <ins>Pipes</ins>


`who`			Display who is logged in

One method to get a sorted list of names is to type,
`who > names.txt`
`sort < names.txt`

This is a bit slow, and you must remember to remove the temporary file called names.txt when you have finished. What you really want to do is connect the output of the who command directly to the input of the sort command. This is exactly what pipes do. The symbol for a pipe is the vertical bar `|`

`who | sort`

- Print the count of sorted list of all users logged in
- Using pipes, display all lines of list1 and list2 containing the letter 'p', and sort the result.
<br>


#### <ins>Creating files</ins>

Files containing texts can be created using text editors such as `nano`, `vi`, `vim` on command line

Create a text file named `draft.txt` in a directory called `thesis` using nano and type the texts as shown on the screenshot below.

<br>
<left><img src="/img/nano-screenshot.png" alt="text editing" width="100%"/></left>
<br>


Once we’re happy with our text, we can press `Ctrl+O` (press
the Ctrl or Control key and, while holding it down, press the O key) to write
our data to disk (we’ll be asked what file we want to save this to:
press Return to accept the suggested default of draft.txt).
<br>


#### <ins>Redirect and save output</ins>

To have the output go to both a file and the screen simultaneously.
`cat science.txt | tee out.stdout.txt`

You can also use `tee` to catch `stderr` as shown in this example:
`rmdir no-dir | tee out.stderr.txt`
<br>


#### <ins>Wildcards</ins>

`*` is a wildcard, which matches zero or more characters. Let’s consider the 
`shell-lesson-data/exercise-data/proteins` directory: `*.pdb` matches 
`ethane.pdb`, `propane.pdb`, and every file that **ends with .pdb**. 

On the other hand, `p*.pdb` only matches `pentane.pdb` and `propane.pdb`, 
because the ‘p’ at the front only matches filenames that **begin with the letter ‘p’**.
<br>



#### <ins>File naming conventions</ins>

A directory is merely a special type of file. So, the rules and conventions for naming files apply also to directories.
In naming files, characters with special meanings such as **/ * & %** , should be
avoided. 
Also, avoid using **spaces** within names. 
The safest way to name a file is to use only **alphanumeric characters**, that is, letters and numbers, together with **_** (underscore) and **.** (dot).


## Examples

<br>
<left><img src="/img/file-naming.png" alt="Naming files" width="100%"/></left>
<br>


#### <ins>File permissions</ins>


In UNIX, there are three types of owners: `user`, `group`, and `others`.
File permissions fall in three categories: `read`, `write`, and `execute`.

<br>
<left><img src="/img/file-permission-01.png" alt="File permissions" width="100%"/></left>
<br>


<br>
<left><img src="/img/file-permission-02.png" alt="File permissions" width="100%"/></left>
<br>


#### <ins>Changing permissions </ins>

`chmod` 					change file modes

**Syntax**					
`chmod permissions filename`
Where,
permissions can be `read`, `write`, `execute` or a combination of them.
filename is the name of the file for which the permissions need to change. 
This parameter can also be a list if files to **change permissions in bulk**.

We can change permissions using **two modes**:
1. **Symbolic mode**: this method uses symbols like `u`, `g`, `o` to represent `users`, `groups`, and `others`. Permissions are represented as  `r`, `w`, `x` for `read`, `write` and `execute`, respectively. You can modify permissions using `+`, `-` and `=`.

2. **Absolute mode**: this method represents permissions as **3-digit octal
   numbers ranging from 0-7**.


#### <ins>Symbolic mode</ins>
Suppose I have a script `fastqc-run.sh` and want to make it executable for the
**owner**, `jjuma`.

<br>
<left><img src="/img/fastqc-run-01.png" alt="File permissions" width="100%"/></left>
<br>

I can change the file permission using the command

`chmod u+x fastqc-run.sh`

<br>
<left><img src="/img/fastqc-run-02.png" alt="File permissions" width="100%"/></left>
<br>

What would the commands below do to the script?

- `chmod go+x fastqc-run.sh`
- `chmod o-r fastqc-run.sh`
- `chmod g=w fastqc-run.sh`
 
#### <ins>Absolute mode</ins>
Absolute mode uses numbers to represent permissions and mathematical operators to modify them.

<br>
<left><img src="/img/absolute-mode-01.png" alt="File permissions" width="100%"/></left>
<br>


1. Set **read** for **user**, **read** and **execute** for **group**, and only **execute** for others.

<br>
<left><img src="/img/absolute-mode-03.png" alt="File permissions" width="100%"/></left>
<br>

`chmod 451 filename`

2. Remove **execution** rights from **others** and **group**.

<br>
<left><img src="/img/absolute-mode-04.png" alt="File permissions" width="100%"/></left>
<br>

`chmod 440 filename`


**Exercise**
If the resulting mode of `chmod 451` is equivalent to `r--r-x--x`?

Write the resulting file modes for the following
- `chmod 777 filename`
- `chmod 754 filename`
- `chmod 700 filename`











