Assignment 7 Due March 20 at 10am

Please provide the exact command line arguments you used to generate your results.
{How to run demultiplex.py}
-
Question 1:
{Control Count}
    - Control: 19
{Treatment Count}
    - Treatment: 24
{Non-Determined Count}
    - Non-determined cells: 7

Question 2:
{Control Count}
    - Control: 19 + 3 = 22
{Treatment Count}
    - Treatment = 24 + 2 = 26
{Non-Determined Count}
    - Non-determined cells: 7 - 5 = 2

Question 3:
{3 strategies to maximize group assignment}
    - Extend the length of the barcode as much as possible without increasing average of errors in tags-with-errors. Then, relax hamming distance threshold
    - Test different locations for the tag -- maybe there are locations less prone to mutation
    - Whitelist tags with hamming distance 1 (or <= the threshold). Then, accept tags that have greater hamming distance
      from known cellTag, but a hamming distance of <= 1 (or some threshold) to that of a cellTag on the whitelist

Comments:
{Things that went wrong or you can not figure out}
-
Suggestions:
{What programming and/or genomics topics should the TAs cover in the next class that would have made this assignment go smoother?}
