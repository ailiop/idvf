## Contributing to idvf


We welcome contributions that help improve the `idvf` project!  Please read
through the guidelines in this document before reporting an issue or
submitting a request.  We will do our best to respond to all issues and
requests, but please bear in mind that it may take us a while.



<a name="toc"></a>

### Contents


- [Types of contributions](#contrib-types)
- [Bug reports](#bug-reports)
- [Pull requests](#pull-requests)
- [Code and documentation style](#style)
- [Feature requests](#feature-requests)



<a name="contrib-types"></a>

### Types of contributions


-   *Bug reports.*
-   *Compatibility patches:* support for earlier MATLAB versions;
    alternatives to MATLAB Toolbox functions.
-   *Minor patches:* documentation clarification and typos;
    code/documentation formatting; naming convention amendments.
-   *Functionality updates.*
-   *Testing:* testing or demo scripts for existing functionality.
-   *Re-implementations:* performance improvements; system/language support
    extensions.
-   Anything else, as long as its utility and functionality is described.



<a name="bug-reports"></a>

### Bug reports


Please [open a new issue][github-new-issue] for each unreported bug.
Specify "[BUG] *1-sentence-description-of-bug*" as the issue title, and
list the following information in the issue body:

-   Brief summary and background.
-   Bug description: what should happen, and what happens instead.
-   Version of MATLAB and relevant toolboxes; operating system; GPU
    hardware and drivers (if the bug affects GPU execution).
-   Code for a concise MATLAB script that reproduces and illustrates the
    bug.
-   Any other relevant notes (e.g., what you think causes the bug, any
    steps you may have taken to identify or resolve it, etc).


[github-new-issue]: https://help.github.com/articles/creating-an-issue/



<a name="pull-requests"></a>

### Pull requests


Please [submit a GitHub pull request][github-pull-request] for each code or
documentation contribution to `idvf`.  When submitting a pull request,
please adhere to the following.

-   Clearly identify the [type of your contribution](#contrib-types) in the
    title and body of your pull request.
    -   If your contributions span multiple types, please separate them
        into individual pull requests.  Minor patches should be lumped into
        a single pull request.
-   Include a brief description of the rationale, functionality, and
    implementation of your contribution.
-   Include a testing or demo MATLAB script (named `test_xxx.m` or
    `demo_xxx.m`) that can be used to illustrate and validate your
    contribution.  Include a brief description of the script in the pull
    request body.
-   [Squash partial commits][github-squash-commit].
-   If applicable, draft some relevant text to be added to or amended in
    the README.  Please include the text in the pull request body, *not* as
    part of the commit.

We encourage you to open a new issue to discuss any intended contributions
prior to developing or submitting a pull request.


[github-pull-request]:  https://help.github.com/articles/about-pull-requests/

[github-squash-commit]: https://help.github.com/articles/about-pull-request-merges/



<a name="style"></a>

### Code and documentation style


Please try to follow the style conventions in the `idvf` repository when
submitting pull requests.  We recommend that you use the source code in
`+dvf/inversion.m` as a style template for functions, and
`demo_inversion_2d.m` for scripts.  We generally try to observe the
following rules:

-   The code should be clear, stable, and efficient.  Clarity and stability
    take precedence over efficiency and performance.  The code should be
    self-documented if possible (avoid referring to descriptions in
    existing issues or pull requests).
-   Function documentation should be comprehensive and follow the format of
    existing functions (e.g., `+dvf/inversion.m`).
-   Function and variable names are in `camelCase`; script names are in
    `snake_case`.  Generally, matrix/array names start with an uppercase
    letter, while scalar/vector/function names start with a lowercase
    letter.
-   All code blocks should be briefly documented.
-   We prefer 4-space indentation (no tabs), operator/operand alignment
    across multiple lines, and 80-column line width.



<a name="feature-requests"></a>

### Feature requests


It is unlikely that we will do much development for new features, unless
they are essential or supported by theoretical advances.  We do, however,
encourage the development and submission of new features via pull requests.
