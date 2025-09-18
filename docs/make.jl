using Documenter, Obs
push!(LOAD_PATH, "../src")
makedocs(sitename = "Obs Documentation", doctest=true,
	 repo = "https://github.com/AntoninoDAnna/Obs.git",
	 pages = [
		      "Home" => "index.md",
		      "Enum" => "enum.md",
		      "Improvements" => "improvements.md",
		      "Observables" =>  "obs.md",
		      "Utilities" => "utilities.md"],
	 format = Documenter.HTML(repolink = "https://github.com/AntoninoDAnna/Obs.gi",
                                  prettyurls=false))



