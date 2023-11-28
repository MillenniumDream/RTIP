//! About the warning and error information when an interrupt occurs at running time.





/// Error message for File reading, creating, opening, and writing.
pub fn error_file(operation: &str, filename: &str) -> String
{
    format!("\n\n\n ERROR: There is some problem in {} the file '{}'. \n\n\n", operation, filename)
}

/// Error message for Directory creating
pub fn error_dir(operation: &str, dir: &str) -> String
{
    format!("\n\n\n ERROR: There is some problem in {} the directory '{}'. Maybe it already exists or you have no permission. \n\n\n", operation, dir)
}





/// Error message for CString::new() 
pub fn error_str_to_cstring(str_name: &str) -> String
{
    format!("\n\n\n ERROR: There is some problem in transforming str '{}' to CString. \n\n\n", str_name)
}

/// Error message for try_into()
pub fn error_type_transformation(variable: &str, type1: &str, type2: &str) -> String
{
    format!("\n\n\n ERROR: There is some problem in type transformation of '{}' from {} to {}. \n\n\n", variable, type1, type2)
}

/// Error message for as_slice() and as_slice_mut()
pub fn error_as_slice(variable: &str) -> String
{
    format!("\n\n\n ERROR: There is some problem in getting the slice of the variable '{}'. \n\n\n", variable)
}





/// Error message for `Some<A>`, Result<T, E>
pub fn error_none_value(variable: &str) -> String
{
    format!("\n\n\n ERROR: There is some problem with variable '{}', which has none/wrong value. \n\n\n", variable)
}

/// Error message for cloned()
pub fn error_cloning(variable: &str) -> String
{
    format!("\n\n\n ERROR: There is some problem in cloning the variable '{}'. \n\n\n", variable)
}





/// Error message for min_1d function
pub fn error_min_1d() -> String
{
    format!("\n\n\n ERROR: There is some problem with the function min_1d: the input fun is increasing along +x direction, or the default minimum step is too large. \n\n\n")
}

/// Error message for illegal chemical element type
pub fn error_element() -> String
{
    format!("\n\n\n ERROR: Illegal chemical element type has read from the input file. Please check it. \n\n\n")
}

/// Error message for read_xyz function
pub fn error_read_xyz(filename: &str) -> String
{
    format!("\n\n\n ERROR: There is some problem with the input file '{}'. Please check it. \n\n\n", filename)
}










