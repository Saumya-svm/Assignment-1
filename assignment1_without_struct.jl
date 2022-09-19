using Printf

const R = 8.314

global temp

global p = 1000

println("**Kindly enter the correct temperature range**")
println("The code below has assumed no temperature constraints to exist. This implies that even though the low temperature propertyx of N2 is 300, the code will run even on temperatures below the low temperature")

"""
# Description:
- Extracting the properties from the first line

# Arguments:
- line from reading the data stream from the therm.dat file

# Return:
- None
"""

function extract_meta_data(line)
    line = split(line, "  ")
    keys = reverse(["low_temp", "high_temp", "common_temp"])
    global species 
    species = string(line[1])
    values = []
    try
        values = (reverse(line)[3:5])
        values = map(a->parse(Float64, a), values)
    catch
        return ""
    end
    species_properties[species] = Dict{Any, Any}(zip(keys, values))
    species_properties[species]["high_temp_properties"] = Dict()
    species_properties[species]["low_temp_properties"] = Dict()
    species_coefficients[species] = Dict("high_temp"=>Dict(), "low_temp"=>Dict())
end

"""
# Description
- Extract coefficients from the second line for upper temperature interval

# Arguments
- line from reading the data stream from the therm.dat file 

# Return
- None

# Exception
- None
"""
function second_line(line)
    i = 1
    while i < 76
        try
            species_coefficients[species]["high_temp"]["a"*string(floor(Int, i/15 +1))] = parse(Float64,line[i:i+14])
        catch
            return ""
        end
        i+=15
    end
end


"""
# Description
- Extract coefficients from the second line for upper and lower temperature interval

# Arguments
- line from reading the data stream from the therm.dat file

# Return
- None

# Exception
- None
"""

function third_line(line)
    i = 76
    while i < 106
        try
            species_coefficients[species]["high_temp"]["a"*string(floor(Int, i/15 +1))] = parse(Float64,line[i-75:i-75+14])
        catch
            return ""
        end
        i+=15
    end
    i = 1
    while i < 46
        try
            species_coefficients[species]["low_temp"]["a"*string(floor(Int, i/15 +1))] = parse(Float64,line[i+30:i+44])
        catch
            return ""
        end
        i+=15
    end

end


"""
# Description
- Extract coefficients from the second line for lower temperature interval

# Arguments
- line from reading the data stream from the therm.dat file

# Return
- None

# Exception
- None
"""
function fourth_line(line)
    i = 46
    while i < 106
        try
            species_coefficients[species]["low_temp"]["a"*string(floor(Int, i/15 +1))] = parse(Float64,line[i-45:i-45+14])
        catch
            return ""
        end
        i+=15
    end
end


"""
# Description
- Extract coefficients according to the temperature

# Arguments
- T(Float64) - Temperature at which coeffs have to be extracted

# Return
- coeff(Dict) - dictionary containing the coeff

# Exception
- 
"""
function get_coefficients(T::Float64, species::String)
    try
        if T > species_properties[species]["common_temp"] 
            coeff = species_coefficients[species]["high_temp"]
            temp_range = "high_temp_properties"
        else
            coeff  = species_coefficients[species]["low_temp"]
            temp_range = "low_temp_properties"
        end
        return coeff
    catch
        println("Enter the correct temperature range")
        T = parse(Float64, readline())
        temp = T
        return get_coefficients(temp, species, true)
    end
end

"""
# Description:
- This function can be used to calculate the speicifc heat of a species.

# Arguments:
- T :: Float64
- species(optional)::String

# returns:
- Specific heat value at the given temperature.

"""
function calculate_specific_heat(T::Float64, species=species) 
    temp = T
    coeff = get_coefficients(temp, species)
    species_specific_heat = R*(coeff["a1"] + coeff["a2"]*temp + coeff["a3"]* (temp^2) + coeff["a4"]*(temp^3) + coeff["a5"]*(temp^4))
    species_properties[species]["specific_heat"] = species_specific_heat
    return species_specific_heat
end

"""
# Description:
- This function can be used to calculate the enthalpy of a species.

# Arguments:
- T :: Float64
- species(optional)::String

# returns:
- Enthalpy value at the given temperature.

"""
function calculate_enthalpy(T::Float64, species=species)
    temp = T
    coeff= get_coefficients(T, species)
    species_enthalpy = R*temp*(coeff["a1"] + coeff["a2"]*(temp)/2 + coeff["a3"]*(temp^2)/3 + coeff["a4"]*(temp^3)/4 + coeff["a5"]*(temp^4)/5 + coeff["a6"]/temp)
    species_properties[species]["enthalpy"] = species_enthalpy
    return species_enthalpy
end

"""
# Description:
- This function can be used to calculate the entropy of a species.

# Arguments:
- T :: Float64
- species(optional)::String

# returns:
- entropy value at the given temperature.

"""
function calculate_entropy(T::Float64, species=species)
    temp = T
    coeff = get_coefficients(T, species)
    species_entropy = R*(coeff["a1"]*log(temp) + coeff["a2"]*(temp) + coeff["a3"]* (temp^2)/2 + coeff["a4"]*(temp^3)/3 + coeff["a5"]*(temp^4)/4 + coeff["a7"])
    species_properties[species]["entropy"] = species_entropy
    return species_entropy
end




# reading a file
f = open("therm.dat", "r")

# creating and empty dictionary to store the species data
species_properties = Dict()
species_coefficients = Dict()

# reading the first line just to ignore the blank line at the top
readline(f)

# there are other ways to parse the data as well. Those include - 
# 1. Reading the lines one by one and only parse those who have a length of around 82 as mentioned in the pdf.
# 2. Assume the data will have a line gap at the top and read the line before looping through the file.
# Here I have used the second method however method one could have been used. The index for the function dictionary could have been accesses through the line number
# at the end of each line. This parsing code is according to the pdf provided. 

lines  = readlines(f)
function_dict = Dict(1=>extract_meta_data, 2=>second_line, 3=>third_line, 0=>fourth_line)
for (i, line) in enumerate(lines)
    function_dict[i%4](line)
end

# till now we have tabulated the coefficients for all the available speices

"""
This function could be used for calculating properties of all species at once.
"""
function calculate_species_properties(T::Float64)
    for i in keys(species_properties)
        global species
        species = i
        species_properties[species]["specific_heat"] = calculate_specific_heat(T)   
        species_properties[species]["enthalpy"] = calculate_enthalpy(T)
        species_properties[species]["entropy"] = calculate_entropy(T)    
    end
end



function calculate_individual_species_property(Species::String, T::Float64)
    global species
    species = Species
    temp = T
    species_properties[species]["specific_heat"] = calculate_specific_heat(temp)
    species_properties[species]["enthalpy"] = calculate_enthalpy(temp)
    species_properties[species]["entropy"] = calculate_entropy(temp)
    println("\nThermodynamic properties of pure $Species at $T(K) and $p(Pa)\n")
    println("Species \t Enthalpy(J/mol) \t Entropy(J/mol-K) \t cp(J/mol-K)\n")
    @printf("%s \t\t %+.4e \t\t %+.4e \t\t %.4e \n", species, species_properties[species]["enthalpy"], species_properties[species]["entropy"], species_properties[species]["specific_heat"])
    # return species_properties[species]["specific_heat"], species_properties[species]["enthalpy"], species_properties[species]["entropy"]
end


"""
# Description
- Calculates the specific heat of a mixture of species

# Arguments
- mixture_comp(Dict) : contains the fractions of each species in the mixture
- T(Float64) : 

# Return
- returns the mixture specific heat

# Exception
- KeyError
"""
function calculate_mixture_specific_heat(mixture_comp, T)
    mixture_specific_heat = 0.0
    global species
    temp = T
    for i in keys(mixture_comp)
        species = i
        species_specific_heat= calculate_specific_heat(temp, species)
        species_properties[species]["specific_heat"] = species_specific_heat
        mixture_specific_heat += species_specific_heat*mixture_comp[species]
    end
    return mixture_specific_heat
end


"""
# Description
- Calculates the enthalpy of a mixture of species

# Arguments
- mixture_comp(Dict) : contains the fractions of each species in the mixture
- T(Float64) : 

# Return
- returns the mixture enthalpy

# Exception
- KeyError
"""
function calculate_mixture_enthalpy(mixture_comp, T)
    mixture_enthalpy = 0.0
    global species
    temp = T
    for i in keys(mixture_comp)
        species = i
        species_enthalpy= calculate_enthalpy(temp, species)
        species_properties[species]["enthalpy"] = species_enthalpy
        mixture_enthalpy += species_enthalpy*mixture_comp[species]
    end
    return mixture_enthalpy
end


"""
# Description
- Calculates the entropy of a mixture of species

# Arguments
- mixture_comp(Dict) : contains the fractions of each species in the mixture
- T(Float64) : 

# Return
- returns the mixture specific heat

# Exception
- KeyError
"""
function calculate_mixture_entropy(mixture_comp, T)
    mixture_entropy = 0.0
    global species
    temp = T
    for i in keys(mixture_comp)
        species = i
        species_entropy= calculate_entropy(temp, species)
        species_properties[species]["entropy"] = species_entropy
        mixture_entropy += species_entropy*mixture_comp[species]
    end
    return mixture_entropy
end


"""
# Description
- Calculate mixture properties Cp, H, S

# Arguments
- mixture_comp(Dict) : contains the fractions of each species in the mixture
- T(Float64) : 

# Return
- 

# Exception
- 
"""
function CalculateMixtureProperties(mixture_comp::Dict, T::Float64)
    temp = T
    mixture_specificheat= calculate_mixture_specific_heat(mixture_comp, temp)
    mixture_enthalpy= calculate_mixture_enthalpy(mixture_comp, temp)
    mixture_enntropy= calculate_mixture_entropy(mixture_comp, temp)

    p = 1000
    println("\nThermodynamic properties of pure species at $T(K) and $p(Pa)\n")
    println("Species \t Enthalpy(J/mol) \t Entropy(J/mol-K) \t cp(J/mol-K)\n")
    l = ["CH4", "CO", "CO2", "H2", "H2O", "O2", "N2"]
    for species in l
        # println(species)
        if mixture_comp[species] == 0.0
            species_properties[species]["specific_heat"] = calculate_specific_heat(temp, species)  
            species_properties[species]["enthalpy"] = calculate_enthalpy(temp, species)
            species_properties[species]["entropy"] = calculate_entropy(temp, species)
            @printf("%s \t\t %+.4e \t\t %+.4e \t\t %.4e \n", species, species_properties[species]["enthalpy"], species_properties[species]["entropy"], species_properties[species]["specific_heat"])
        else
            @printf("%s \t\t %+.4e \t\t %+.4e \t\t %.4e \n", species, species_properties[species]["enthalpy"], species_properties[species]["entropy"], species_properties[species]["specific_heat"])
        end
    end
    println()
    println("The mixture enthalpy is ", mixture_enthalpy)
    println("The mixture entropy is ", mixture_enntropy)
    println("The mixture specific heat is ", mixture_specificheat)
end