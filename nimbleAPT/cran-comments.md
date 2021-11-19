## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a new release.

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

The example I provide for buildAPT takes 30 seconds to run. A key bottle neck is that NIMBLE models and sampling schemes are typically compiled prior to running MCMC.
So I have replaced \dontrun with \donttest.


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
