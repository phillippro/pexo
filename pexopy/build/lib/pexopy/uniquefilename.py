import hashlib

class UniqueFilename(str):
   """
   Generates an MD5 hash from the `contents` (<str>) and returns a string in a format f'{prepend}{md5}{append}'
   """
   def __new__(cls, contents, prepend="", append="", *args, **kwargs):
      md5 = hashlib.md5(contents.encode("utf-8")).hexdigest()
      name = f"{prepend}{md5}{append}"
      # TODO : handle existing files?
      return str.__new__(cls, name)