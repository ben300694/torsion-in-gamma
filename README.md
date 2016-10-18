# Torsion in (Gamma pi_2 K) / pi_1 K

This code was written for the
bachelor thesis in mathematics with the topic
"Torsion in ![equation](http://www.sciweavers.org/tex2img.php?eq=%5CGamma%20%28%5Cpi_2%20K%29%20%2F%20%5Cpi_1%20K&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)".

The PDF can be found [here](https://github.com/ben300694/torsion-in-gamma/blob/master/pdfs/20160605_bachelorthesis_torsion_in_gamma.pdf).

Author:
Benjamin Ruppik (University of Bonn)

Date:
June 8th, 2016

Supervisors:
Dr. Daniel Kasprowski and
Prof. Dr. Peter Teichner
(Max Planck Institute for Mathematics in Bonn)

This is the Git-repository for publishing the
bachelor thesis and the SageMath code.

## Instructions

You will need [SageMath](http://www.sagemath.org/index.html) to run this program (version 7.2 or above is recommended).

The relevant SageMath code is contained in the directory `src/`.
Navigate there in your terminal, this directory should contain the file `main.sage`.

Now start SageMath in command line mode by typing `sage`,
then load the main file via `attach('main.sage')`.

The input for most of the relevant functions is the presentation of a finite, finitely presented group.
See the corresponding [SageMath documentation](http://doc.sagemath.org/html/en/reference/groups/sage/groups/finitely_presented.html)
on how to enter such a presentation.

The function you probably want to use is
```python
test(grouppresentation, grouppresentation_string="default", progress=True, use_c_code=False)
```

Suppose `G` is a finitely presented group, then you would call this function as
```python
test(G, "meaningful_name_for_G")
```
The function call creates two new files in the directory:
- `meaningful_name_for_G_output.txt` contains the results.
- `meaningful_name_for_G_progress.txt` is updated live to reflect the progress the program has made.

For example if you want to test the standard presentation of the cyclic group with 5 elements
you would call
```python
test(cyclic_group(5), "Z5")
```

For more examples see the bachelor thesis.

**The program might take a very long time to process some group presentations, please be patient.**
You can use the terminal command `tail -f "meaningful_name_for_G_progress.txt"`
to check up on the progress of the calculation.

## Plans for the future

- [x] Implement basic multithreading functionality
- [ ] Write critical parts of the program in C or another more efficient language to speed up the process
