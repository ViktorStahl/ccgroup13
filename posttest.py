import json 
def postResult(result):
	feeds = {}
	try:
		with open('result.txt', mode='r') as file:
			feeds={}
			try:
				feeds = json.load(file)
			except:
				pass
	except:
		pass
	with open('result.txt', mode='w+') as file:
		feeds.update(result)
		json.dump(feeds, file)
		file.close()
		return '0'

if __name__ == '__main__':
	postResult({'1COS': {'problem': '1', 'time': '2.41', 'method': 'COS', 'result': '0.923'}})