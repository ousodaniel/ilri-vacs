# Linux for Bioinformatics

---
***Trainer: OUSO | [LinkedIn](https://github.com/ousodaniel) | [GitHub](https://github.com/ousodaniel)***
## Course Overview

This course assumes familiarity with the Linux operating system (OS). It brushes over the basics, while progressively advancing. It focuses on skills required for bioinformatics analysis.

By the end of this course, participants will be:
- Develop basic scripts for automation
- Work from high-performance computing clusters using SLURM

### Target Audience
- Basic Linux users
- Researchers, students, and professionals in bioinformatics and computational biology
- Anyone who needs to use Linux-based tools for analysing biological data

### Prerequisites
- Familiarity with Linux operating system
- Interest in bioinformatics applications

## Linux Fundamentals

---
### Introduction to Linux
---
**Brief History and Philosophy**

Linux was created in 1991 by Linus Torvalds as a free, open-source alternative to Unix. The Linux operating system is built around the Linux kernel, with various software components added to create complete distributions. The philosophy behind Linux emphasises:

- Open-source development
- Community collaboration
- Freedom to use, modify, and distribute
- Stability and security
- Transparency

**Linux Distributions (Distros)**

Linux comes in many forms called _distributions_ or _distros_. Common ones include:

- **Ubuntu**: User-friendly, good for beginners
- **CentOS/RHEL**: Common on enterprise servers and HPC clusters
- **Debian**: Known for stability
- **Fedora**: Cutting-edge features
- **Scientific Linux/Biolinux**: Specialized for scientific computing

For bioinformatics, **Ubuntu** and **CentOS** are widely used, with many specialised bioinformatics tools pre-configured to work well on these platforms.

**Key Differences from Windows**

| Feature | Windows | Linux |
|---------|---------|-------|
| User Interface | Primarily GUI-based | Command-line focused (with GUI options) |
| File System | Drive letters (`C:`, `D:`) | Unified hierarchical structure (`/`) |
| Case Sensitivity | Not case-sensitive | Case-sensitive |
| File Extensions | Required | Optional |
| Terminal/Command Line | PowerShell/CMD | Bash, Zsh, etc. |
| Software Installation | `.exe` installers | Package managers |
| System Customization | Limited | Highly customizable |
| Cost | Licensed | Free (mostly) |

**Why Linux for Bioinformatics?**

Linux dominates bioinformatics for several reasons:

1. **Performance**: Efficient with computational resources
2. **Command Line Interface**: Enables automation and reproducible analyses
3. **Tool Availability**: Most bioinformatics tools are developed for Linux first
4. **Scripting Capabilities**: Easy to create pipelines and workflows
5. **Remote Access**: Simple `SSH` connectivity to high-performance computing resources
6. **Open Source Nature**: Aligns with scientific principles of transparency and reproducibility

#### Accessing Linux from Windows

**Windows Subsystem for Linux (WSL)**

WSL lets you run a Linux environment directly on Windows without a standalone virtual machine.

*Setting up WSL:*
1. Open PowerShell as Administrator
2. Run: `wsl --install`
3. Restart your computer
4. Launch Ubuntu from the Start menu
5. Create a username and password

<details>
    <summary>
        Click to toggle contents of <b style='color:blue'>Virtual Machines</b>
    </summary>
    A virtual machine (VM) creates a Linux environment within your Windows computer.

*Setting up VirtualBox with Ubuntu:*
1. Download VirtualBox from virtualbox.org
2. Download an Ubuntu ISO file
3. Create a new VM in VirtualBox
4. Allocate RAM (4GB minimum recommended for bioinformatics)
5. Create a virtual hard disk (50GB+ recommended)
6. Start the VM and select the ISO file to install Ubuntu
</details>

<details>
    <summary>Click to toggle contents of <b style='color:blue'>SSH to Remote Servers</b>
    </summary>
Many bioinformatics analyses run on remote Linux servers or clusters, usually HPCs.

*Connecting via SSH:*
1. In Windows, install PuTTY or use OpenSSH in PowerShell
2. For OpenSSH, use: `ssh username@server-address`
3. Enter your password when prompted
</details>


#### The Linux File System

A file system (FS) is the method and structure an operating system uses to store, organise, and manage data on storage devices.

**File System Hierarchy**

Linux organises files in a tree-like structure starting from the root directory[^1] (`/`):
[^1]: _Folder_ equivalent in Windows

<details>
    <summary>
        Click to toggle contents of <b style='color:blue'>Full Linux file structure </b>
    </summary>
May not be updated.
    
```
/
├── bin/    # Essential commands
├── boot/   # Boot loader files
├── dev/    # Device files
├── etc/    # System configuration
├── home/   # User home directories
├── lib/    # Shared libraries
├── media/  # Removable media
├── mnt/    # Temporary mounts
├── opt/    # Optional software
├── proc/   # Process information
├── root/   # Root user's home
├── run/    # Run-time data
├── sbin/   # System administration commands
├── srv/    # Service data
├── tmp/    # Temporary files
├── usr/    # User utilities and applications
└── var/    # Variable files (logs, etc.)
```
    
</details>

**Important Directories**

- `/home`: Where user files are stored (similar to `C:\Users` in Windows)
- `/usr`: Contains programs, libraries, and documentation
- `/etc`: System-wide configuration files
- `/bin`: Essential commands
- `/sbin`: System binaries
- `/tmp`: Temporary files (cleared on reboot)
- `/var`: Data that changes during system operation (logs, databases)

**Path Concepts**

- **Absolute paths** start from the root directory: `/home/username/data/sequences.fa`
- **Relative paths** start from the current directory: `data/sequences.fa`
- **Home directory** shorthand: `~` (tilde)
- **Parent directory**: `..` (double fullstop)
- **Current directory**: `.` (single fullstop)
- **Previous directory**: `-` (hyphen)

---
### Command Line Basics
---
#### Concepts

**Kernel:** The core operating system component running in the background. It manages hardware resources; handles memory, CPU, and device management and allows communication between software and hardware.
 
**Shell:** The user interface for accessing the services of the kernel. It interprets the command line inputs and calls the appropriate programs or system calls.

**Terminal:** The interface that allows users to interact with the shell. It can be a graphical application like GNOME Terminal or a virtual console (`tty`).

#### Terminal Essentials

**Opening the Terminal**
- In Ubuntu: `Ctrl+Alt+T` or search for `Terminal`
- In WSL: Open `Ubuntu` from the `Start menu`

**Command Syntax**

Commands generally take the below form:
```
command [options] [arguments]
```

However, sometimes it is possible to get things done with the `command` alone, or with the `command` and either of `options` or `arguments`. The square brackets around the latter two components indicate that they are optional.
- `command`: The program to run (usually the name of the program in lowercase with hyphens for spaces)
- `options`: Modify/Specify the command's behaviour (usually preceded by `-` (short form) or `--` (long form))
- `arguments`: What the command acts on: files, directories, links etc

**Using Tab Completion**

`Tab` completion is a shell feature (used in Linux, macOS, and other Unix-like systems) that auto-completes commands, file names, directory paths, or other inputs when you press the `Tab` key on your keyboard.
- Start typing a command or file name and press `Tab`
- The shell will complete it if there's only one possibility
- Press `Tab` twice to see all possibilities

**Command History**

We write code because we are lazy and prefer not to unnecessarily repeat tasks. Linux stores your command line history in a register to allow quick retireval for repeat use cases--saved time.
- `Up/Down` arrows to navigate through previous commands
- `history` to list recent commands: note the register number of your command and rerun it by preceding the number with a `!` (exclamation).
- `Ctrl+R` to search command history

**Getting Help**

Help is never far from you, at least one of the following options should work for your case:
- `man <command>`: Display the manual page for a command
- `<command> -h`: Brief help for a command
- `<command> --help`: Extended help for a command
- `info <command>`: Detailed documentation (if available)

#### Basic Commands
Basic Linux commands are inspired by natural language (English), thus may be easy to work with.

**Navigation Commands**
- `pwd`: Print working directory (shows current location)
- `ls`: List files and directories
  - `ls -l`: Detailed listing
<!--   - `ls -a`: Show hidden files
  - `ls -h`: Human-readable file sizes -->
- `cd <directory_name>`: Change directory
    - `cd -`: Go to previous directory
<!--   - `cd` or `cd ~` - Go to home directory
  - `cd ..`: Go up one directory (parent dir) -->
  

**Directory and File Operations**
- `mkdir <directory_name>`: Create directory
  - `mkdir -p parent/child`: Create parent directories as needed
- `touch <file>`: Create empty file or update timestamp
- `cp <source> <destination>`: Copy files or directories
  - `cp -r <source> <destination>`: Copy directories recursively
- `mv <source> <destination>`: Move or rename files/directories
- `rm <file_name>`: Remove files
  - `rm -r <directory_name>`: Remove directories recursively
  - `rmdir <directory_name>`: Remove empty directory
<!--   - `rm -f file`: Force removal without confirmation -->


**Viewing File Contents**
- `cat <file_name>`: Display entire file
- `less <file_name>`: View file with pagination (use `q` to quit)
- `head <file_name>`: Show first 10 lines
  - `head -n 20 <file_name>`: Show first 20 lines
- `tail <file_name>`: Show last 10 lines
  - `tail -n 20 <file_name>`: Show last 20 lines
  - `tail -f <file_name>`: Follow file updates in real-time

**System Information**
- `whoami`: Display current username
- `date`: Show current date and time
- `echo <text>`: Display text
- `uname -a`: Show system information

---
### Working with Files and Directories
---

#### File Manipulation

**Text Editors**

**Nano** (beginner-friendly):
- Open/Create: `nano <file_name>`
- Save: `Ctrl+O`, then `Enter`
- Exit: `Ctrl+X`
- Cut line: `Ctrl+K`
- Paste: `Ctrl+U`
- Find text: `Ctrl+W`

**File Permissions**

Linux uses a permission system for files and directories:

```
-rwxrwxrwx
 ↑  ↑  ↑
 |  |  └─ Others (everyone else)
 |  └──── Group (users in the file's group)
 └─────── Owner (user who owns the file)
```

Each set has three permission types:
- `r` (read): View file contents or list directory contents
- `w` (write): Modify file or create/delete files in directory
- `x` (execute): Run file as program or access files in directory

**Changing Permissions:**
`chmod <permissions> <file>`
<!-- - Numeric method: read=4, write=2, execute=1
  - `chmod 755 script.sh` (rwxr-xr-x)
  - `chmod 644 data.txt` (rw-r--r--) -->
- Symbolic method:
  - `chmod u+x script.sh` (to user [`u`] add [`+`] execute [`x`] rights)
  - `chmod go-w data.txt` (to group [`g`] and others [`o`] remove [`-`] write [`w`] rights)

<!-- **Changing Ownership:**
- `chown user:group file`
- Example: `chown john:researchers data.fa` -->

<!-- **Hidden Files and Directories**
- Names start with a dot (.)
- Not shown by default with `ls`
- View with `ls -a`
- Common examples: `.bashrc`, `.profile`, `.cache` -->

#### Finding and Organising Files

**The `find` Command**

The `find` command locates files based on various criteria:

```bash
# Find all .fastq files in the current directory and subdirectories
find . -name "*.fastq"

# Find and delete temporary files
find ./temp -name "*.tmp" -delete

# Find and execute a command on each file
find ./data -name "*.txt" -exec cat {} \;
```

**The `grep` command for Searching Within Files**

`grep` searches for patterns within files line-wise:

```bash
# Search for lines with "ATGC" in a file
grep "ATGC" sequence.fa

# Search for lines without "ATGC" in a file
grep -v "ATGC" sequence.fa

# Count occurrences
grep -c "ATGC" sequence.fa

# Show only file_names with matches
grep -l "ATGC" *.fa

# Show context (3 lines before and after)
grep -B 3 -A 3 "ATGC" sequence.fa
```

**Wildcards and Pattern Matching**

- `*`: Matches any number of characters
  - `*.fastq`: All files with `.fastq` extension
  - `sample_*`: All files starting with `sample_`
- `?`: Matches a single character
  - `sample?.txt`: Matches `sample1.txt`, `sampleA.txt`, etc.
- `[]`: Matches any character within brackets
  - `sample[123].txt`: Matches `sample1.txt`, `sample2.txt`, or `sample3.txt`
- `{}`: Group of patterns
  - `{*.fastq,*.fq}`: Matches files ending in `.fastq` or `.fq`

---
### Linux Redirection and Pipes
---

#### Input/Output Redirection

What happens when you want to feed the output of one command as the inout to another? What about if you want to save the output of a command, whether that output is an _error_ or the _results_ from a process?

**Standard Streams**

Linux commands use three standard streams:
- **`stdin (0)`**: Standard input (keyboard by default)
- **`stdout (1)`**: Standard output (terminal by default)
- **`stderr (2)`**: Standard error (terminal by default)

**Redirection Operators**

- `>`: Redirect `stdout` to a file (***overwrites***)
  ```bash
  # Save directory listing to a file
  ls -l > file_list.txt
  ```

- `>>`: Redirect `stdout` to a file (***appends***)
  ```bash
  # Append current date to a log file
  date >> log.txt
  ```

- `2>`: Redirect `stderr` to a file
  ```bash
  # Save error messages to a file
  find / -name "*.fa" 2> errors.log
  ```

- `2>&1`: Redirect `stderr` to same destination as `stdout`
  ```bash
  # Save both output and errors
  ls /bin /nonexistent > output.txt 2>&1
  ```

- `<`: Use file as `stdin`
  ```bash
  # Use file as input
  sort < unsorted.txt
  ```

<!-- - `<<`: Provide multi-line input (document) realtime; notice End Of File marker at the end of input
  ```bash
  cat << EOF > sequences.fa
  >Sequence1
  ATGCTAGCTAGCTAGCT
  >Sequence2
  GCATGCATGCATGCAT
  EOF
  ``` -->

#### Pipes, Filters, and Preprocessing

**The Pipe Operator**

The pipe (`|`) sends the output of one command as input to another:

```bash
# List files and filter for FASTQs; files ending ($) with .fastq
ls | grep "\.fastq$"

# Count number of FASTA files; files ending ($) with .fa
ls | grep "\.fa$" | wc -l

# Count number of sequences in a FASTA file; lines starting (^) with > (greater than)
grep "^>" sequences.fa | wc -l
```

**Common Filter Commands**

- `sort`: Sort lines of text
  ```bash
  # Sort a file alphabetically
  sort gene_names.txt
  
  # Sort numerically
  sort -n measurements.txt
  
  # Sort by column (tab-delimited)
  sort -k2,2n data_table.tsv
  ```

- `uniq`: Report or filter out repeated lines
  ```bash
  # Find unique lines (input must be sorted)
  sort gene_names.txt | uniq
  
  # Count occurrences of each line
  sort gene_names.txt | uniq -c
  ``` 

- `wc`: Count lines, words, and characters
  ```bash
  # Count lines
  wc -l file.txt
  
  # Count words
  wc -w file.txt
  
  # Count characters
  wc -c file.txt
  ```

- `cut`: Extract columns from files
  ```bash
  # Extract 2nd column from tab-delimited file
  cut -f2 data.tsv
  
  # Extract 1st and 3rd columns from CSV
  cut -d',' -f1,3 data.csv

**Preprocessing**

- `tr`: Character-level Replace, Delete, Collapse
  ```bash
  # Replace
  tr "," "\t" < data.csv > data.tsv
  
  # Delete
  tr -d "\r" < data-carriage-returns.csv > data.csv
  
  # Collapse
  echo "gene     chromosome     position" | tr -s ' ' 
  ```
  
- `sed`: String-level processing
  ```bash
  # Substitution
  sed 's/ATGC/NNNN/' <<< "ATGC ATGC"
  
  sed 's/ATGC/NNNN/g' <<< "ATGC ATGC"
  
  # Print lines with matching <pattern>
  sed -n '/pattern/p' <file.txt>
  
  # Delete
  sed '/^$/d' <file.txt>
  sed '10,20d' file.txt # delete lines 10 - 20
  
  #Anchoring
  sed 's/^>/HEADER:/' # prefix
  
  sed 's/$/_processed/' # suffix
  ```
  
- `awk`: Handling structured data; line-by-line and column-wise
  ```bash
  # awk 'CONDITION/PATTERN { ACTION }'
  # Print columns
  awk -F'\t' '{ print $1, $2 }' data.tsv
  
  # Mean of a column
  awk '{ sum += $2 } END { print sum / NR }' data.tsv
  
  # BEGIN and END
  # Run before processing input
  awk 'BEGIN { print "Starting..." } { print $1 }' data.tsv
  
  # Run after processing all input
  awk '{ count++ } END { print count }' data.tsv
  
  # Use Bash variables
  min_length=500

  awk -v min="$min_length" '$2 >= min' data.tsv
  ```
<details>
    <summary>
        Click to toggle contents of <b style='color:blue'>AWK Quick Reference </b>
    </summary>
May not be updated.
    
| Expression  | Meaning                 |
| ----------- | ----------------------- |
| `$0`        | Entire line             |
| `$1`        | First field             |
| `$NF`       | Last field              |
| `NF`        | Number of fields        |
| `NR`        | Current line number     |
| `/pattern/` | Lines matching regex    |
| `==`        | String/numeric equality |
| `~`         | Regex match             |
| `!~`        | Regex non-match         |
| `-F'\t'`    | Tab delimiter           |
| `BEGIN {}`  | Before input            |
| `END {}`    | After input             |

</details>

---
### Shell Scripting Basics
---
#### Introduction to Shell Scripts
A script is a file containing a sequence of commands for task automation.

**Creating and Executing Scripts**

1. Create a script file: `my_script.sh`
   ```bash
   #!/bin/bash
   # My first bioinformatics script
   
   echo "Starting analysis..."
   echo "Current directory: $(pwd)"
   echo "Files to process:"
   ls -1 *.fa
   ```

2. Make it executable (if you don't provide any of the `u`, `g` and `o` permission qualifiers, the right is granted to everyone):
   ```bash
   chmod +x my_script.sh
   ```

3. Run it:
   ```bash
   ./my_script.sh
   ```

**Shebang Line**

The first line of a script (`#!/bin/bash`) tells the system which interpreter to use:
- `#!/bin/bash`: Bash shell
- `#!/bin/sh`: POSIX shell (more portable)
<!-- - `#!/usr/bin/env python`: Python
- `#!/usr/bin/env Rscript`: R script -->

**Variables and Environment**

***Variable***: a named symbol that stores value which can be evaluated in a command or script.

***Environment***: a collection of key-value pairs that define the behaviour of processes.

```bash
# Set variables
project_name="Genome Analysis"
data_dir="/path/to/data"
num_samples=10

# Use variables
echo "Project: $project_name"
echo "Analyzing $num_samples samples from $data_dir"

# Command substitution
today=$(date +%Y-%m-%d)
echo "Report date: $today"

# Environment variables
echo "User home: $HOME"
echo "Username: $USER"
echo "Path: $PATH"

# Reading user input
echo -n "Enter sequence file: "
read file_name
echo "You selected: $file_name"
```

#### Script Control Structures

**Conditionals**

What if you want to do something based on whether some logic evaluates to `true` or `false`?
```bash
#!/bin/bash

# If-else statement
in_file="sequence.fa"
if [ -f "$in_file" ]; then
    echo "$in_file exists."
else
    echo "$in_file does not exist."
    exit 1
fi

# Check file size
file_sz=$(stat -c%s "$in_file")
if [ $file_sz -gt 1000000 ]; then
    echo "$in_file is large (> 1MB)."
elif [ $file_sz -gt 1000 ]; then
    echo "$in_file is medium (> 1KB)."
else
    echo "$in_file is small."
fi

# Test command options
# -f: regular file exists
# -d: directory exists
# -z: string is empty
# -n: string is not empty
# -eq, -ne, -lt, -le, -gt, -ge: numeric comparison
```

**Loops**

How can we perform the same operations over some sequence of input?
```bash
#!/bin/bash

# For loop with files: exhaustive
for file in *.fastq; do
    echo "Found FASTQ: $file"
    file_bs=$(basename "$file" .fastq)
    echo "Base name: $file_bs"
done

# While loop: trigger-based
count=1
while [ $count -le 5 ]; do
    echo "Count: $count"
    count=$((count+1))
done
```

**Command Substitution**

It is possible to substitute the output of a command in another command.
```bash
# Store command output in variable
num_seqs=$(grep -c "^>" sequence.fa)
echo "Found $num_seqs sequences"

# Use in expressions
if [ $(wc -l < data.txt) -gt 100 ]; then
    echo "File has more than 100 lines"
fi
```

#### Bioinformatics Script Example

```bash
#!/bin/bash
# Script: extract_long_seqs.sh
# Purpose: Extract sequences longer than specified length

# Check arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 <fasta_file> <min_length>"
    exit 1
fi

fasta_file=$1
min_length=$2

# Validate input file
if [ ! -f "$fasta_file" ]; then
    echo "Error: File '$fasta_file' not found."
    exit 1
fi

# Create output file_name
out_file="${fasta_file%.fa}_min${min_length}.fa"
echo "Extracting sequences longer than $min_length bp from $fasta_file"
echo "Output will be saved to $out_file"

# Process the file
# Initialize variables
current_header=""
current_seq=""
count=0
total=0

# Process each line
while IFS= read -r line; do
    # If header line
    if [[ $line == ">"* ]]; then
        # Process previous sequence if we have one
        if [ -n "$current_seq" ]; then
            total=$((total+1))
            length=${#current_seq}
            if [ $length -ge $min_length ]; then
                echo "$current_header" >> "$out_file"
                echo "$current_seq" >> "$out_file"
                count=$((count+1))
            fi
        fi
        # Start new sequence
        current_header=$line
        current_seq=""
    else
        # Append to current sequence
        current_seq="${current_seq}${line}"
    fi
done < "$fasta_file"

# Process the last sequence
if [ -n "$current_seq" ]; then
    total=$((total+1))
    length=${#current_seq}
    if [ $length -ge $min_length ]; then
        echo "$current_header" >> "$out_file"
        echo "$current_seq" >> "$out_file"
        count=$((count+1))
    fi
fi

echo "Done. Extracted $count out of $total sequences."
```

**Running the Script:**

```bash
chmod +x extract_long_seqs.sh
./extract_long_seqs.sh sequences.fa 1000
```

---
### Introduction to High-Performance Computing with SLURM
---

#### HPC Concepts

**High-Performance Computing (HPC)**

High-Performance Computing refers to using computing clusters for computationally intensive tasks. Key benefits for bioinformatics include:

- Processing large datasets that exceed desktop capabilities
- Running analyses in parallel
- Access to more memory, storage, and computing power
- Specialized hardware for specific tasks

**Cluster Architecture**

A typical HPC cluster consists of:

- **Login node(s)**: Where users connect and manage jobs
- **Compute node(s)**: Where jobs actually run
- **Storage systems**: For data and results
- **Interconnect**: High-speed network fabric connecting nodes
- **Job scheduler**: Software that manages job execution

**Resource Management Concepts**

- **Jobs**: Computational tasks submitted to the cluster
- **Resources**: CPUs, memory, time, GPUs, etc.
- **Queues/Partitions**: Job categories with different priorities and resources
- **Allocation**: Amount of resources available to a user/group
- **Walltime**: Maximum execution time for a job

#### SLURM Basics

**What is SLURM?**

SLURM (Simple Linux Utility for Resource Management) is a workload manager for Linux clusters. It provides three key functions:

1. Allocate resources to users
2. Provide a framework for starting/monitoring jobs
3. Manage a queue of pending jobs

**Key SLURM Commands**

- `sbatch`: Submit a batch job script
  ```bash
  sbatch my_job.sh
  ```

- `squeue`: View information about jobs in the queue
  ```bash
  # Show all jobs
  squeue
  
  # Show only your jobs
  squeue -u $USER
  ```

- `scancel`: Cancel a job
  ```bash
  # Cancel a specific job
  scancel 12345
  
  # Cancel all your jobs
  scancel -u $USER
  ```

- `sinfo`: View information about compute nodes and partitions
  ```bash
  # Show cluster status
  sinfo
  
  # Show detailed node information
  sinfo -N -l
  ```

#### Package management
THe HPC  is a shared space. User requirements may vary substantially, there is a centralised way to ensure module/software/package needs are met from a central place (eliminate/minimise redundancy). Availability of modules is managed using the `module` command.

```bash
# List available module
module list

# Load modules for use
module load <module_name> [<module_name1> <module_name2> ...]

# Which modules are available for your use
module avail

# Unload modules
module unload <module_name> [<module_name1> <module_name2> ...]

# Reset module availability
module purge
```

#### Creating SLURM Job Scripts

**Basic SLURM Script Structure**

```bash
#!/bin/bash
#SBATCH --job-name=myanalysis      # Job name
#SBATCH --output=myanalysis_%j.out # Output file (%j expands to job ID)
#SBATCH --error=myanalysis_%j.err  # Error file
#SBATCH --time=01:00:00            # Wall time (HH:MM:SS)
#SBATCH --nodes=1                  # Number of nodes
#SBATCH --ntasks=1                 # Number of tasks
#SBATCH --cpus-per-task=4          # Number of CPUs per task
#SBATCH --mem=8G                   # Memory per node
#SBATCH --partition=normal         # Partition/queue name

# Load any required modules
module load samtools bwa

# Run commands
echo "Job started at $(date)"
echo "Running on host: $(hostname)"
echo "Working directory: $(pwd)"

# Example bioinformatics commands
bwa mem -t $SLURM_CPUS_PER_TASK reference.fa read1.fq read2.fq > output.sam
samtools sort -@ $SLURM_CPUS_PER_TASK -o output.bam output.sam
samtools index output.bam

echo "Job finished at $(date)"
```

<details>
    <summary>
        Click to toggle contents of 
        <b style='color: blue'> SLURM Directives
        </b>
    </summary>
    Common SLURM directives for bioinformatics jobs:

| Directive | Description | Example |
|-----------|-------------|---------|
| `--job-name` | Name of the job | `--job-name=alignment` |
| `--output` | Standard output file | `--output=align_%j.out` |
| `--error` | Standard error file | `--error=align_%j.err` |
| `--time` | Maximum wall time | `--time=12:00:00` |
| `--nodes` | Number of nodes | `--nodes=1` |
| `--ntasks` | Number of tasks | `--ntasks=1` |
| `--cpus-per-task` | CPUs per task | `--cpus-per-task=8` |
| `--mem` | Memory per node | `--mem=32G` |
| `--mem-per-cpu` | Memory per CPU | `--mem-per-cpu=4G` |
| `--partition` | Partition to use | `--partition=compute` |
| `--account` | Account to charge | `--account=bioinformatics` |
| `--mail-user` | Email for notifications | `--mail-user=user@example.com` |
| `--mail-type` | When to send email | `--mail-type=END,FAIL` |
| `--array` | Job array indices | `--array=1-10` |
| `--dependency` | Job dependencies | `--dependency=afterok:12345` |

</details>

**Environment Variables**

SLURM sets environment variables that can be used in your script:

```bash
echo "Job ID: $SLURM_JOB_ID"
echo "Number of nodes: $SLURM_JOB_NUM_NODES"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Node list: $SLURM_JOB_NODELIST"
```

---
### Bioinformatics Workflows with SLURM
---

#### Running Bioinformatics Tools with SLURM

**Example: DNA Sequence Alignment**

```bash
#!/bin/bash
#SBATCH --job-name=bwa_align
#SBATCH --output=bwa_align_%j.log
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

# Load modules
module load bwa samtools

# Set variables
ref_genome="/path/to/reference/genome.fa"
sample_id="sample1"
r1="/path/to/reads/${sample_id}_R1.fastq.gz"
r2="/path/to/reads/${sample_id}_R2.fastq.gz"
out_dir="/path/to/output"

# Create output directory
mkdir -p $out_dir

# Start timing
start_time=$(date +%s)

# Run alignment
echo "Starting alignment for $sample_id"
bwa mem -t $SLURM_CPUS_PER_TASK $ref_genome $r1 $r2 | \
    samtools sort -@ $SLURM_CPUS_PER_TASK -o $out_dir/${sample_id}.bam -

# Index BAM file
samtools index $out_dir/${sample_id}.bam

# Calculate statistics
samtools flagstat $out_dir/${sample_id}.bam > $out_dir/${sample_id}.flagstat

# End timing
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Job completed in $elapsed_time seconds"
```

## Additional Resources

### Online Documentation and Tutorials

- **Linux Command Line**:
  - [The Linux Command Line](http://linuxcommand.org/) - comprehensive guide
  - [Software Carpentry Shell Lessons](https://swcarpentry.github.io/shell-novice/)
  - [Explain Shell](https://explainshell.com/) - interactive command explanation

- **Bioinformatics-specific Resources**:
  - [Bioinformatics Workbook](https://bioinformaticsworkbook.org/)
  - [Biostars Handbook](https://www.biostarhandbook.com/)
  - [UNIX for Biologists](https://unix4biologists.org/)

- **SLURM Documentation**:
  - [Official SLURM Documentation](https://slurm.schedmd.com/documentation.html)
  - [SLURM Cheatsheet](https://slurm.schedmd.com/pdfs/summary.pdf)

## Glossary of Terms

| Term | Definition |
|------|------------|
| **Absolute Path** | A file path that starts from the root directory (`/`) |
| **Argument** | Information passed to a command |
| **Bash** | Bourne Again SHell, the most common Linux shell |
| **CPU** | Central Processing Unit, performs computations |
| **Directory** | A folder in the file system |
| **Environment Variable** | A named value that can affect how processes run |
| **File System** | The structure used to control how data is stored and retrieved |
| **Flag** | A type of option that modifies command behavior |
| **GUI** | Graphical User Interface |
| **HPC** | High-Performance Computing |
| **Job** | A computational task submitted to a cluster |
| **Kernel** | The core of the operating system |
| **Node** | A single computer in a cluster |
| **Option** | A modifier for a command |
| **Path** | The location of a file or directory |
| **Pipe** | A mechanism to pass output from one command to another |
| **Process** | A running program |
| **Redirection** | Changing where input comes from or output goes to |
| **Relative Path** | A file path that starts from the current directory |
| **Shell** | A command-line interpreter |
| **SLURM** | Simple Linux Utility for Resource Management |
| **Standard Error** | The default stream for error messages |
| **Standard Input** | The default stream for input |
| **Standard Output** | The default stream for program output |
| **Terminal** | A program that runs a shell |
| **Variable** | A named storage location |
| **Wildcard** | A character that represents multiple characters |

