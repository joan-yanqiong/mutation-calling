output_dir="/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/output/test_set_2"
mkdir -p $output_dir
while read symlink_path actual_path
do
    # mkdir -p $(dirname $symlink_path)
    # mv $symlink_path $(dirname $actual_path)
    # echo $actual_path
    cp -r $actual_path $output_dir
done < "/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/output/tumor_true_file_paths.txt"