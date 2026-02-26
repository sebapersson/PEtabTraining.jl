using PEtabTraining
using Documenter
using DocumenterVitepress

DocMeta.setdocmeta!(PEtabTraining, :DocTestSetup, :(using PEtabTraining); recursive = true)

format = Documenter.HTML(;
    prettyurls = get(ENV, "CI", "false") == "true",
    assets = String["assets/custom_theme.css"],
    repolink = "https://github.com/sebapersson/PEtabTraining.jl",
    edit_link = "main"
)

makedocs(;
    modules = [PEtabTraining],
    sitename = "PEtabTraining.jl",
    repo = Remotes.GitHub("sebapersson", "PEtabTraining.jl"),
    authors = "Sebastian Persson, and contributors",
    checkdocs = :exports,
    warnonly = false,
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/sebapersson/PEtabTraining.jl",
    ),
    pages = [
        "Home" => "index.md",
        "API" => "API.md",
        "Contributing" => "CONTRIBUTING.md",
    ],
)

DocumenterVitepress.deploydocs(
    repo = "github.com/sebapersson/PEtabTraining.jl.git",
    target = "build", # this is where Vitepress stores its output
    devbranch = "main",
    branch = "gh-pages",
    push_preview = true,
)
