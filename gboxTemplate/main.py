import granatum_sdk


def main():
    gn = granatum_sdk.Granatum()

    assay = gn.get_import('assay')
    assay_df = gn.pandas_from_assay(assay)

    checkbox_value = gn.get_arg('someCheckbox')
    number_value = gn.get_arg('someNumber')
    seed_value = gn.get_arg('someSeed')

    markdown_str = f"""\
  * checkbox_value = **{checkbox_value}**
  * number_value = **{number_value}**
  * seed_value = **{seed_value}**
  * Shape of the assay = **{assay_df.shape}**"""

    gn.add_result(markdown_str, data_type='markdown')

    gn.commit()


if __name__ == "__main__":
    main()
