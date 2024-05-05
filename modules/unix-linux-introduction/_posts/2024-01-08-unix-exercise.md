---
title: Hands-on exercise
---

## Hands-on exercise in UNIX
<br>

1. Download the `cases.csv` file using `wget` command

```
wget https://hpc.ilri.cgiar.org/~jjuma/training_data/cases.csv
```


The `cases.csv` file contains some simulated data on cholera cases in Haiti. Once eradicated
from the island, cholera reappeared in the Department (or region) of Artibonite in 2010 and
more than 800,000 cases of disease have been recorded over the past decade.


If you've got this far, you've found that you've been provided with a file on cholera case data, from
the island state of Haiti. 

{% figure [caption:"Departments of Haiti"] [class:"caption"] %}
![Nine Departments of Haiti with their names: Artibonite, Centre, Grande'Anse, Nippes, Nord, Nord-Est, Nord-Ouest, Sud, Sud-Est, Ouest](/img/Haiti_departments_named.png)
{% endfigure %}

1. How many cases are recorded in the file? 
<details markdown="1">
<summary>Hint</summary>
To answer this question, you need to extract the departments column (the 2nd column) and find out how many unique values are in that column.
</details>
2. What are the headings? 
3. In the case ID, why was the location where the case was collected not encoded as part of the ID?


<details markdown="1">
There are two common ways to extract a column from a file, one using `cut` and one using `awk`. In both cases these tools are designed with pipes in mind, so
they can either read from a file or from `stdin`. And in both cases you need to take note of what the field separator is: in this case it is a comma (","),
making this file a Comma Separated Values (CSV). Another common file type that you might encounter is a TSV, where columns are separated with tabs.

Of the two tools, `cut` is the simpler, and it uses `-d` to set the delimiter (which is a tab by default) and `-f` to select a field:

```bash
cut -d, -f2 cases.tsv
```

The `awk` tools is more powerful but also more complex. It is very well worth learning how to use `awk` and there are good guides online (including some of those
listed as part of our course resources). For `awk` we also need to specific a field separator, with `-F`, and then a small program to extract the field we want.
The program `{ print $2}` prints the 2nd field for each line in the file. Notice that it uses variables with `$` signs, just like the shell, so to stop
the shell interpreting our `awk` program we have to specific it in single quotes:

```bash
awk -F, '{print $2}'
```

If you look closely, you'll notice that both `cut` and `awk` shows us the first line of the file (the header) and we don't want to see that. While we could use some
`awk` programming to avoid that, the simplest way to start reading from the second line onwards is using `tail` with the `-n+2` flag. If you specify locations with
a `+`, `tail` will start reading from that line onwards, instead of counting backwards from the end of the file. A pipe lets us combine the two commands:

```bash
tail -n+2 cases.csv |cut -d, -f2
```
<summary>Show answer</summary>
</details>

The first part of the answer shows us how to extract the data we need from the file. Now, how to to count the unique values?

<details markdown="1">

In the shell, `sort` and `uniq` let us turn a list of values into a list of distinct values. And `wc` lets us count, so these three programs form part
of the everyday vocabulary of shell users.

The `uniq` command only works on repeated lines, i.e. lines that are right next to each other. This is why we need to ensure that data is sorted with `sort`
before using `uniq`. 

<summary>Show answer</summary>

Putting it all together:

```bash
tail -n+2 cases.csv |cut -d, -f2 |sort |uniq
```

will show us the unique values in column 2 and:

```bash
tail -n+2 cases.csv |cut -d, -f2 |sort |uniq|wc -l
```

will count how many unique departments there are in our dataset (there are 9).
</details>

The third column of the dataset is the date that the case was collected. This is organised in YEAR-MONTH-DAY order with 0-padded months and days, which
aligns with the [ISO 8601](https://en.wikipedia.org/wiki/ISO_8601) standard. One advantage of this way of writing dates is that they are easily sorted.

Using this information, how can we find the maximum date recorded in the file?

<details markdown="1">
<summary>Show answer</summary>
The tool to use here is `sort`. On its own `sort` can focus on a particular column (known as a `key`). Once again we need to tell the tool how to separate
columns, this time using the `-t` flag. So as a first attempt:

```bash
sort -t, -k3 cases.csv
```

works to sort the data by date, from earliest to latest. Again, we probably want to use `tail -n+2` to skip the header.

There is another thing to note about `-k3`: the default sort order used by `sort` is dictionary sort. For this data, that doesn't make a difference, but
if we want to sort numbers, we need to use `-k3n`. If our dates were expressed without the leading zeros it would cause a problem for sort. Luckily our
ISO-8601-compliant dates don't have this problem.

Again we might want to combine a few commands:

```bash
tail -n+2 cases.csv |cut -d, -f3 |sort |tail -1
```

What answer does that give us?
</details>

To find the minimum date, we need a few small changes? What are they

<details markdown="1">
<summary>Show answer</summary>

We could either reverse the sort order with `-r`:

```bash
tail -n+2 cases.csv |cut -d, -f3 |sort -r |tail -1
```

or we could pick the first line from the original smallest-to-largest sort:

```bash
tail -n+2 cases.csv |cut -d, -f3 |sort |head -1
```

</details>

Shifting focus again, how can we find out which department has the most cases? 

<details markdown="1">
<summary>Show answer</summary>
The `uniq` tool has a `-c` flag for counting how often distinct values appear. We can use that here:

```bash
tail -n+2 cases.csv |cut -d, -f2 | sort | uniq -c
```

Remember, always `sort` before `uniq`!
</details>

Male and female gender is encoded as M or F in the dataset. Thank goodness that this is clean data! Which gender saw the most cases?

<details markdown="1">
<summary>Show answer</summary>
We can use a very similar set of commands to answer this question, now just using column 4:

```bash
tail -n+2 cases.csv |cut -d, -f4 | sort | uniq -c
```
</details>

Finally, how many cases were seen in children under 1 year of age? We have to look carefully at the data to answer this.

<details markdown="1">
<summary>Show answer</summary>
The age is expressed in months for children under age 1. We can use that to filter out the data we need using the `grep` command:

```bash
tail -n+2 cases.csv |cut -d, -f6 | grep months | wc -l
```

Or we could use `awk`. Once again `awk` is a bit more complex but more flexible:

```bash
tail -n+2 cases.csv |awk -F, '$6 ~ /months/ {tot = tot + 1} END {print tot}'
```

The test used in `awk` here is actually a regular expression (i.e. pattern) match so we can use it to check for more than just the word "months".
</details>

## Resources

[The Carpentries Shell Lesson](https://swcarpentry.github.io/shell-novice/)

[Terminal Basics from Sandbox.Bio](https://sandbox.bio/tutorials/terminal-basics)

[Bash Scripting](https://www.linode.com/docs/guides/intro-bash-shell-scripting/)

[Data Exploration with Awk](https://sandbox.bio/tutorials/awk-intro)

[The Baeldung Linux Guides](https://www.baeldung.com/linux/)

[Grep, Sed and Awk](https://www.linode.com/docs/guides/differences-between-grep-sed-awk/)

[Bashcrawl - an adventure in Bash exploration](https://gitlab.com/slackermedia/bashcrawl)

[Linux Cheat Sheet](https://www.geeksforgeeks.org/linux-commands-cheat-sheet/)


