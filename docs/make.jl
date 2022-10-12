using Documenter
using Apolo

makedocs(
    sitename = "Apolo",
    format = Documenter.HTML(),
    modules = [Apolo]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/mvanzulli/Apolo.jl.git"
)
