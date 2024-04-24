## Changes in version 1.0.5 - 23/04/2024
- A bug fix for the frequency of printed output during the adaptation of the temperature ladder.
- Removed conflict in printed output between temperature ladder and progress bar.


## Changes in version 1.0.5 - 22/04/2024
A simple change, added resetMV as a logical flag in the run function of buildAPT.
This addition makes the package more in line with NIMBLE.


## Results of 'R CMD check --as-cran nimbleAPT_1.0.4.tar.gz'
* using log directory ‘/home/pleydell/nimbleProjects/nimbleAPT/nimbleAPT.Rcheck’
* using R version 4.1.2 (2021-11-01)
* using platform: x86_64-pc-linux-gnu (64-bit)
* using session charset: UTF-8
* using option ‘--as-cran’
* checking for file ‘nimbleAPT/DESCRIPTION’ ... OK
* checking extension type ... Package
* this is package ‘nimbleAPT’ version ‘1.0.4’
* package encoding: UTF-8
* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘David Pleydell <david.pleydell@inrae.fr>’

New submission
* checking package namespace information ... OK
* checking package dependencies ... OK
* checking if this is a source package ... OK
* checking if there is a namespace ... OK
* checking for executable files ... OK
* checking for hidden files and directories ... OK
* checking for portable file names ... OK
* checking for sufficient/correct file permissions ... OK
* checking serialization versions ... OK
* checking whether package ‘nimbleAPT’ can be installed ... OK
* checking installed package size ... OK
* checking package directory ... OK
* checking for future file timestamps ... OK
* checking ‘build’ directory ... OK
* checking DESCRIPTION meta-information ... OK
* checking top-level files ... OK
* checking for left-over files ... OK
* checking index information ... OK
* checking package subdirectories ... OK
* checking R files for non-ASCII characters ... OK
* checking R files for syntax errors ... OK
* checking whether the package can be loaded ... OK
* checking whether the package can be loaded with stated dependencies ... OK
* checking whether the package can be unloaded cleanly ... OK
* checking whether the namespace can be loaded with stated dependencies ... OK
* checking whether the namespace can be unloaded cleanly ... OK
* checking loading without being on the library search path ... OK
* checking use of S3 registration ... OK
* checking dependencies in R code ... OK
* checking S3 generic/method consistency ... OK
* checking replacement functions ... OK
* checking foreign function calls ... OK
* checking R code for possible problems ... OK
* checking Rd files ... OK
* checking Rd metadata ... OK
* checking Rd line widths ... OK
* checking Rd cross-references ... OK
* checking for missing documentation entries ... OK
* checking for code/documentation mismatches ... OK
* checking Rd \usage sections ... OK
* checking Rd contents ... OK
* checking for unstated dependencies in examples ... OK
* checking installed files from ‘inst/doc’ ... OK
* checking files in ‘vignettes’ ... OK
* checking examples ... OK
* checking for unstated dependencies in vignettes ... OK
* checking package vignettes in ‘inst/doc’ ... OK
* checking re-building of vignette outputs ... OK
* checking PDF version of manual ... OK
* checking for non-standard things in the check directory ... OK
* checking for detritus in the temp directory ... OK
* DONE

Status: 1 NOTE
See
  ‘/home/pleydell/nimbleProjects/nimbleAPT/nimbleAPT.Rcheck/00check.log’
for details.


__So just one note saying that this is a new submission.__


## In response to Gregor Seyer email, "CRAN submission nimbleAPT 1.0.4", 18 November 2021 15:10

> Please add \value to .Rd files regarding exported methods and explain
> the functions results in the documentation. Please write about the
> structure of the output (class) and also what the output means. (If a
> function does not return a value, please document that too, e.g.
> \value{No return value, called for side effects} or similar)
> Missing Rd-tags:
>       buildAPT.Rd: \value
>       samplers.Rd: \value

Done.

> \dontrun{} should only be used if the example really cannot be executed
> (e.g. because of missing additional software, missing API keys, ...) by
> the user. That's why wrapping examples in \dontrun{} adds the comment
> ("# Not run:") as a warning for the user.
> Does not seem necessary.
>
> Please unwrap the examples if they are executable in < 5 sec, or create
> additionally small toy examples to allow automatic testing.
> (You could also replace \dontrun{} with \donttest, if it takes longer
> than 5 sec to be executed, but it would be preferable to have automatic
> checks for functions. Otherwise, you can also write some tests.)

The example I provide for buildAPT was taking 30 seconds to run.
A key bottle neck is that NIMBLE models and sampling schemes are typically compiled prior to running MCMC.
So I have removed all compilation steps from the examples and refer the user to the vignette for more details.


> You write information messages to the console that cannot be easily suppressed.
> It is more R like to generate objects that can be used to extract the
> information a user is interested in, and then print() that object.
> Instead of print()/cat() rather use message()/warning()  or
> if(verbose)cat(..) (or maybe stop()) if you really have to write text to the console.
> (except for print, summary, interactive functions)

In the 'setup' code of buidAPT I have repaced 'print' with 'messages'.
In the 'run' code of buildAPT this is not possible, because run code must be compiled using the NIMBLE compiler.
So here nimPrint must be used instead of R functions such as message().
All such messages in the 'run' code have flags that can be used to turn off this printing.
I have removed two non-essential printed messages from the run code.


> Please make sure that you do not change the user's options, par or
> working directory. If you really have to do so within functions, please
> ensure with an *immediate* call of on.exit() that the settings are reset
> when the function is exited. e.g.:
> ...
> oldpar <- par(no.readonly = TRUE)    # code line i
> on.exit(par(oldpar))            # code line i + 1
> ...
> par(mfrow=c(2,2))            # somewhere after
> ...
> e.g.: APT_functions.R
> If you're not familiar with the function, please check ?on.exit. This
> function makes it possible to restore options before exiting a function
> even if the function breaks. Therefore it needs to be called immediately
> after the option change within a function.

Done.


## In response to Gregor Seyer email, "CRAN submission nimbleAPT 1.0.3", 11 September 2020 05:28

> Please omit the space within the doi specification to make it clickable._

This has been done.

> Please add small executable examples in your Rd-files to illustrate the
> use of the exported function but also enable automatic testing._

More examples have been added.


> Please always add all authors, contributors and copyright holders in the
> Authors@R field with the appropriate roles.
>  From CRAN policies you agreed to:
> "... Where code is copied (or derived) from the work of others
> (including from R itself), care must be taken that any copyright/license
> statements are preserved and authorship is not misrepresented.
> Preferably, an ‘Authors@R’ would be used with ‘ctb’ roles for the
> authors of such code. Alternatively, the ‘Author’ field should list
> these authors as contributors.
> Where copyrights are held by an entity other than the package authors,
> this should preferably be indicated via ‘cph’ roles in the ‘Authors@R’
> field, or using a ‘Copyright’ field (if necessary referring to an
> inst/COPYRIGHTS file)."
> e.g.: Daniel Turek
> Please explain in the submission comments what you did about this issue._

I have contacted the NIMBLE Development Team regarding this issue and they have listed four authors who contibuted to the original NIMBLE code which I modified to make this package.
So I have added an Authors@R section in the DESCRIPTION file.
This now lists myself (as the package developer) and the four members of the NIMBLE development team.
These contributions have been cited using "cph". Whilst Gregor Seyer of CRAN suggest "ctb", the NIMBLE development team have suggested (via personal communication) that "cph" would be an appropriate alternative since they do not consider themselves authors of nimbleAPT.
I have included a inst/COPYRIGHTS file which cites the appropriate lines from the inst/COPYRIGHTS file of the NIMBLE pacakge.
In NIMBLE, only C++ code is licenced under GPL (>=2), therefore, since I have only adapted R code (and not C++ code) from NIMBLE I only need to cite the BSD_3_clause licence in inst/COPYRIGHTS.
I have also added a NOTE to my DESCRIPTION file that clarifies the contributions of each contributor.
I hope you find these modifications accepatable.
