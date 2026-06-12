# Function to clean the terminal
function Clean_Terminal()

    # If system == Windows
    if Sys.iswindows()
        Base.run(`cmd /c cls`)

    # If system is based on Unix
    else
        Base.run(`clear`)
    end
    
end