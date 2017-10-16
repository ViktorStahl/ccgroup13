import json 
def postResult(result):
	with open('/home/ubuntu/result.txt', mode='a+') as file:
		feeds={}
		try:
			feeds = json.load(file)
		except:
			pass
		feeds.update(result)
		json.dump(feeds, file)
		file.close()
		return '0'

if __name__ == '__main__':
	postResult({'1COS': {'problem': '1', 'time': '2.41', 'method': 'COS', 'result': '0.923'}})