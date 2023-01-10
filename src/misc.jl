function log_print(verbose, texts...)
    if (verbose == false) || (verbose == 0)
        return nothing
    elseif verbose == true
        verbose = 1
    end
    println(texts[min(length(texts), verbose)])
end
