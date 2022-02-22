using DirectTranscription
using Documenter

DocMeta.setdocmeta!(DirectTranscription, :DocTestSetup, :(using DirectTranscription); recursive=true)

makedocs(;
    modules=[DirectTranscription],
    authors="Grant Hecht",
    repo="https://github.com/GrantHecht/DirectTranscription.jl/blob/{commit}{path}#{line}",
    sitename="DirectTranscription.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://GrantHecht.github.io/DirectTranscription.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/GrantHecht/DirectTranscription.jl",
    devbranch="main",
)
