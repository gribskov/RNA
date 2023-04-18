import yaml


def yaml_read(filename):
    """---------------------------------------------------------------------------------------------
    read the workflow description from the yaml file

    :param filename: string     yaml workflow description
    :return: dict               workflow as a python dictionary
    ---------------------------------------------------------------------------------------------"""
    fp = open(filename, 'r')
    workflow = yaml.load(fp, Loader=yaml.FullLoader)
    fp.close()
    return workflow


if __name__ == '__main__':
    workflow = yaml_read('workflow.yaml')

    exit(0)
