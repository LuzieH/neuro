using Revise
using neuro

function exportall(mod)
    for n in names(mod, all = true)
        if Base.isidentifier(n) && n ∉ (Symbol(mod), :eval)
            @eval mod export $n
        end
    end
end

exportall(neuro)